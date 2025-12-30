#include "llvm/ADT/SmallVector.h"
#include "llvm/Analysis/CFG.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/LoopPass.h"
#include "llvm/Analysis/PostDominators.h"
#include "llvm/Analysis/ScalarEvolution.h"
#include "llvm/Analysis/ScalarEvolutionExpressions.h"
#include "llvm/Analysis/ValueTracking.h"
#include "llvm/Demangle/Demangle.h"
#include "llvm/IR/Analysis.h"
#include "llvm/IR/BasicBlock.h"
#include "llvm/IR/CFG.h"
#include "llvm/IR/Constants.h"
#include "llvm/IR/Dominators.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/Intrinsics.h"
#include "llvm/IR/PassManager.h"
#include "llvm/IR/Use.h"
#include "llvm/IR/Value.h"
#include "llvm/Pass.h"
#include "llvm/Passes/PassBuilder.h"
#include "llvm/Passes/PassPlugin.h"
#include "llvm/Support/Casting.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Transforms/Utils/BasicBlockUtils.h"
#include "llvm/Transforms/Utils/Cloning.h"
#include <cstddef>
#include <map>
#include <string>

using namespace llvm;

// Command-line flag (works when running opt and your plugin registers options)
static cl::opt<bool> PrefetchIMADebug(
    "prefetch-ima-debug",
    cl::desc("Enable debug output for PrefetchIMA slicing/classification"),
    cl::init(false));

static cl::opt<int>
    PrefetchDistance("ima-prefetch-distance",
                     cl::desc("Prefetch distance (iterations ahead)"),
                     cl::init(8)); // Lower default for reordered graphs

static cl::opt<int> DegreeThreshold(
    "ima-degree-threshold", cl::desc("Minimum loop iterations for prefetching"),
    cl::init(4)); // Lower threshold for more aggressive prefetching

static cl::opt<int>
    BatchSize("ima-batch-size",
              cl::desc("Number of addresses to prefetch per iteration (1-8)"),
              cl::init(1)); // Default 1 to avoid out-of-bounds; set higher for
                            // reordered graphs

static cl::opt<int> LocalityHint(
    "ima-locality",
    cl::desc("Cache locality hint (0=none, 1=low, 2=moderate, 3=high/L1)"),
    cl::init(3)); // Default to L1 for reordered graphs

static cl::opt<std::string> FunctionPattern(
    "ima-function-pattern",
    cl::desc("Function name pattern to match (empty = all with IMA)"),
    cl::init("")); // Empty means match all functions with IMA

static cl::opt<bool> EnableCostModel(
    "ima-enable-cost-model",
    cl::desc("Enable cost model to skip unprofitable prefetches"),
    cl::init(true)); // Enabled by default

static cl::opt<int>
    CostThreshold("ima-cost-threshold",
                  cl::desc("Minimum benefit score (0-100) to insert prefetch"),
                  cl::init(20));

static cl::opt<bool> AdaptiveDistance(
    "ima-adaptive-distance",
    cl::desc("Compute prefetch distance based on loop complexity"),
    cl::init(true));

static cl::opt<int>
    MemoryLatencyCycles("ima-memory-latency",
                        cl::desc("Estimated memory latency in CPU cycles"),
                        cl::init(200));

namespace {

/// \brief Compute estimated benefit of prefetching for a loop.
///
/// Returns a score from 0-100 based on:
/// - Trip count estimate (higher = better)
/// - Number of IMA loads (more loads = more benefit)
/// - Loop depth (deeper = more iterations overall)
///
/// \param L The loop to analyze
/// \param IMALoadCount Number of IMA loads in the loop
/// \param SE Scalar evolution for trip count analysis
/// \returns Benefit score 0-100
int computePrefetchBenefit(Loop *L, int IMALoadCount, ScalarEvolution &SE) {
  int score = 0;

  // Factor 1: Estimated trip count (0-40 points)
  unsigned TripCount = SE.getSmallConstantTripCount(L);
  if (TripCount == 0) {
    // Unknown trip count - assume moderate benefit
    score += 20;
  } else if (TripCount > 1000) {
    score += 40;
  } else if (TripCount > 100) {
    score += 30;
  } else if (TripCount > 10) {
    score += 15;
  } else {
    score += 5; // Very small loop - minimal benefit
  }

  // Factor 2: Number of IMA loads (0-30 points)
  if (IMALoadCount >= 4) {
    score += 30;
  } else if (IMALoadCount >= 2) {
    score += 20;
  } else {
    score += 10;
  }

  // Factor 3: Loop depth (0-30 points)
  unsigned Depth = L->getLoopDepth();
  if (Depth >= 3) {
    score += 30;
  } else if (Depth == 2) {
    score += 20;
  } else {
    score += 10;
  }

  return std::min(score, 100);
}

int estimateLoopCycles(Loop *L) {
  int cycles = 0;

  for (BasicBlock *BB : L->blocks()) {
    for (Instruction &I : *BB) {
      if (isa<LoadInst>(&I)) {
        cycles += 4; // L1 cache hit
      } else if (isa<StoreInst>(&I)) {
        cycles += 3;
      } else if (isa<CallInst>(&I) || isa<InvokeInst>(&I)) {
        cycles += 10;
      } else if (isa<BinaryOperator>(&I)) {
        if (I.getType()->isFloatingPointTy()) {
          cycles += 3; // FP operations
        } else {
          cycles += 1; // Integer operations
        }
      } else if (isa<GetElementPtrInst>(&I)) {
        cycles += 1;
      } else if (isa<ICmpInst>(&I) || isa<FCmpInst>(&I)) {
        cycles += 1;
      } else if (isa<BranchInst>(&I)) {
        cycles += 1;
      } else if (isa<PHINode>(&I)) {
        cycles += 0; // PHI nodes are resolved at compile time
      } else {
        cycles += 1; // Default for other instructions
      }
    }
  }

  return std::max(cycles, 1);
}

int computeAdaptiveDistance(Loop *L, int memoryLatency) {
  int loopCycles = estimateLoopCycles(L);

  // Distance = memory_latency / loop_cycles
  // Clamp to reasonable range [2, 64]
  int distance = memoryLatency / loopCycles;
  distance = std::max(2, std::min(64, distance));

  return distance;
}

static bool isVectorSubscriptCall(const CallInst *CI) {
  if (!CI)
    return false;

  const Function *F = CI->getCalledFunction();
  if (!F)
    return false;

  StringRef MangledName = F->getName();
  if (MangledName.empty())
    return false;

  // Fast mangling pattern check for libstdc++ and libc++
  //   _ZNSt6vector...ixE...
  //   _ZNSt3__16vector...ixE...   (libc++)
  if ((MangledName.contains("St6vector") ||
       MangledName.contains("St3__16vector")) &&
      MangledName.contains("ixE"))
    return true;

  // Fallback: demangle for robustness
  std::string Demangled = llvm::demangle(MangledName.str());
  if ((Demangled.find("std::vector") != std::string::npos ||
       Demangled.find("std::__1::vector") != std::string::npos) &&
      Demangled.find("operator[") != std::string::npos)
    return true;

  return false;
}

class SlicingAnalysis {
  Function &Fn;
  FunctionAnalysisManager &Fam;
  LoopInfo &LInfo;
  ScalarEvolution &Se;

  // Derived debug flag (command-line OR env var)
  bool DebugEnabled;

#define DBG(x)                                                                 \
  do {                                                                         \
    if (DebugEnabled)                                                          \
      errs() << x;                                                             \
  } while (0)

public:
  std::map<Value *, SmallPtrSet<const Instruction *, 32>> Slices;

  /// Constructor: Initializes the slicing analysis for a given function.
  /// Performs per-load analysis, computes backward slices, classifies
  /// irregular memory accesses, and attaches metadata for prefetch candidates.
  SlicingAnalysis(Function &F, FunctionAnalysisManager &AM, LoopInfo &LI,
                  ScalarEvolution &Scev)
      : Fn(F), Fam(AM), LInfo(LI), Se(Scev) {
    // Enable debug if either CLI flag set or environment variable set.
    const char *Env = std::getenv("PREFETCH_IMA_DEBUG");
    DebugEnabled = PrefetchIMADebug || (Env && std::string(Env) == "1");

    DBG("\n==============================\n");
    DBG("[PrefetchIMA] Analyzing Function: " << F.getName() << "\n");
    DBG("==============================\n");

    unsigned LoadCount = 0;

    // Iterate over all instructions in the function.
    for (Instruction &I : instructions(F)) {
      // Only handle Load instructions (potential memory access candidates).
      if (auto *LI = dyn_cast<LoadInst>(&I)) {
        ++LoadCount;
        DBG("\n[Load #" << LoadCount << "] " << *LI << "\n");

        // The pointer being loaded from.
        Value *Addr = LI->getPointerOperand();

        // Compute the static backward slice for this address.
        // Includes both data and control dependencies.
        SmallPtrSet<const Instruction *, 32> Slice =
            computeStaticSlice(Addr, /*IncludeControlDep=*/true);
        Slices[Addr] = Slice;

        // If the address is derived from a GetElementPtr instruction,
        // print its structure.
        if (auto *G = dyn_cast<GetElementPtrInst>(Addr)) {
          DBG("  [GEP] Base pointer: " << *G->getPointerOperand() << "\n");
          for (unsigned I = 1; I < G->getNumOperands(); ++I) {
            DBG("    [GEP idx " << I << "] " << *G->getOperand(I) << "\n");
          }
        }

        // Print slice summary.
        DBG("  [Slice size=" << Slice.size() << "]\n");

        // If debug mode is enabled, dump the full slice to stderr.
        if (DebugEnabled) {
          for (const Instruction *S : Slice)
            errs() << "    " << *S << "\n";
        }

        // Classify this memory access as irregular (IMA), regular, or other.
        StringRef Class = classifyIrregularAccess(LI, Slice, LInfo, Se);

        DBG("  ==> Classified as: " << Class << "\n");

        // If identified as an IMA, attach metadata for later prefetch
        // insertion.
        if (Class == "IMA") {
          LLVMContext &C = F.getContext();
          MDNode *N = MDNode::get(C, MDString::get(C, "IMA"));
          if (!LI->getMetadata("prefetch_ima"))
            LI->setMetadata("prefetch_ima", N);
          DBG("  [Tag] Attached metadata: prefetch_ima\n");
        }
      }
    }

    // Summary message after processing all loads.
    DBG("\n[PrefetchIMA] Finished analyzing "
        << LoadCount << " loads in function " << F.getName() << "\n");
    DBG("==============================\n");
  }

  /// \brief Computes a conservative backward static slice for a given SSA
  /// value.
  ///
  /// Traverses data, control, and structural dependencies starting from `Start`
  /// (typically a memory address operand of a load).
  /// Returns all instructions that may influence the computation of `Start`.
  ///
  /// \param Start - The value to slice backward from.
  /// \param IncludeControlDep - If true, include control dependencies using
  ///                            post-dominator analysis.
  /// \returns A set of instructions forming the static backward slice.
  SmallPtrSet<const Instruction *, 32>
  computeStaticSlice(Value *Start, bool IncludeControlDep = true) {
    // Container for all instructions included in the slice.
    SmallPtrSet<const Instruction *, 32> Slice;
    // Worklist for recursive backward traversal.
    SmallVector<const Value *, 32> Worklist;
    Worklist.push_back(Start);

    // Post-dominator tree (used to infer control dependencies).
    PostDominatorTree *PDT = nullptr;
    if (IncludeControlDep)
      PDT = &Fam.getResult<PostDominatorTreeAnalysis>(Fn);

    // Process worklist until all dependent values are exhausted.
    while (!Worklist.empty()) {
      const Value *V = Worklist.pop_back_val();

      // Only instructions contribute to the static slice.
      if (const Instruction *I = dyn_cast<Instruction>(V)) {
        // Skip already-visited instructions (avoid infinite recursion).
        if (!Slice.insert(I).second)
          continue;

        // === Data dependencies ===
        // Traverse all operands of the instruction (backward use-def chain).
        for (const Use &U : I->operands())
          Worklist.push_back(U.get());

        // === PHI nodes ===
        // Include all incoming values from predecessor blocks.
        if (const PHINode *PN = dyn_cast<PHINode>(I))
          for (unsigned I = 0, E = PN->getNumIncomingValues(); I != E; ++I)
            Worklist.push_back(PN->getIncomingValue(I));

        // === Function calls ===
        // Include all argument values that flow into the call site.
        if (const CallInst *CI = dyn_cast<CallInst>(I))
          for (const Value *Arg : CI->args())
            Worklist.push_back(Arg);

        // === Select instructions ===
        // Add both true and false value branches.
        if (const SelectInst *SI = dyn_cast<SelectInst>(I)) {
          Worklist.push_back(SI->getTrueValue());
          Worklist.push_back(SI->getFalseValue());
        }

        // === GetElementPtr instructions ===
        // Include base pointer and all index operands.
        if (const GetElementPtrInst *G = dyn_cast<GetElementPtrInst>(I)) {
          Worklist.push_back(G->getPointerOperand());
          for (unsigned I = 1, E = G->getNumOperands(); I < E; ++I)
            Worklist.push_back(G->getOperand(I));
        }

        // === Alloca instruction (stack variable) ===
        // Add stores that write to this alloca and their stored values.
        if (const AllocaInst *AI = dyn_cast<AllocaInst>(I)) {
          for (const User *U : AI->users()) {
            if (const StoreInst *SI = dyn_cast<StoreInst>(U)) {
              Worklist.push_back(SI);
              Worklist.push_back(SI->getValueOperand());
            }
          }
        }

        // === Load-from-alloca pattern ===
        // Trace back to stores that define the loaded alloca.
        if (const LoadInst *LI = dyn_cast<LoadInst>(I)) {
          if (auto *Alloca = dyn_cast<AllocaInst>(LI->getPointerOperand())) {
            for (const User *U : Alloca->users()) {
              if (const StoreInst *SI = dyn_cast<StoreInst>(U)) {
                const Value *Stored = SI->getValueOperand();
                Worklist.push_back(Stored);
              }
            }
          }
        }

        // === Control dependencies (optional) ===
        // Approximate control flow influences using the PostDominatorTree.
        if (IncludeControlDep && PDT) {
          const BasicBlock *BB = I->getParent();
          for (const BasicBlock *Pred : predecessors(BB)) {
            const Instruction *Term = Pred->getTerminator();
            if (!Term)
              continue;
            if (!PDT->dominates(BB, Pred))
              Worklist.push_back(Term);
          }
        }
      }
    }

    // Return the complete set of influencing instructions.
    return Slice;
  }

  /// \brief Detects whether a given LoadInst represents an indirect array
  /// access.
  ///
  /// An "indirect array access" (IMA) is a pattern like `A[B[i]]`, where
  /// the index itself is loaded from memory rather than being a simple affine
  /// or loop variable.
  ///
  /// This function checks whether the load's address operand follows this
  /// pattern.
  ///
  /// \param L - The load instruction to analyze.
  /// \returns true if the load uses an index computed indirectly (e.g., through
  /// another load).
  bool isIndirectArrayAccess(const LoadInst *L) const {
    // Sanity check: null load instruction.
    if (!L)
      return false;

    // Skip loads that themselves produce pointers (i.e., not actual memory
    // reads).
    if (L->getType()->isPointerTy())
      return false;

    // Expect the pointer operand to come from a GetElementPtr (GEP)
    // instruction. This indicates array-like addressing (e.g., base + index
    // computation).
    const auto *GEP = dyn_cast<GetElementPtrInst>(L->getPointerOperand());
    if (!GEP)
      return false;

    // Extract the candidate index operand from the GEP.
    // Typically, this is an integer operand that represents an array subscript.
    Value *Index = nullptr;
    for (unsigned I = 1; I < GEP->getNumOperands(); ++I)
      if (GEP->getOperand(I)->getType()->isIntegerTy())
        Index = GEP->getOperand(I);
    if (!Index)
      return false;

    DBG("    [IMA-check] Candidate index: " << *Index << "\n");

    // Remove any pointer casts or integer extension/truncation chains.
    // This normalizes the index value before further analysis.
    Index = Index->stripPointerCasts();
    if (auto *SE = dyn_cast<SExtInst>(Index))
      Index = SE->getOperand(0);
    if (auto *ZE = dyn_cast<ZExtInst>(Index))
      Index = ZE->getOperand(0);
    if (auto *TR = dyn_cast<TruncInst>(Index))
      Index = TR->getOperand(0);

    // Case 1: Index itself is loaded from memory (e.g., i = B[j]; A[i]).
    // This is a strong signal of an indirect access pattern.
    const auto *IdxLoad = dyn_cast<LoadInst>(Index);
    if (IdxLoad && IdxLoad->getType()->isIntegerTy()) {
      DBG("    [IMA-check] Index comes from integer load.\n");
      return true;
    }

    // Case 2 (fallback): Index stored in an alloca variable and defined through
    // a store chain. Example:
    //   store (load B[j]), %tmp
    //   ...
    //   load A[%tmp]
    // This pattern still represents an indirect index.
    if (auto *Alloca = dyn_cast<AllocaInst>(Index)) {
      for (User *U : Alloca->users()) {
        if (StoreInst *SI = dyn_cast<StoreInst>(U)) {
          if (isa<LoadInst>(SI->getValueOperand())) {
            DBG("    [IMA-check] Index derived from store chain via load.\n");
            return true;
          }
        }
      }
    }

    // If neither direct nor alloca-based indirect indexing is detected,
    // the load is not an indirect array access.
    return false;
  }

  /// \brief Checks whether a computed static slice contains an index-dependent
  /// load.
  ///
  /// This function inspects the backward slice of a load address (as computed
  /// by `computeStaticSlice`) to determine whether it depends on another memory
  /// load that produces an integer value.
  ///
  /// A load of an integer value within the slice indicates *index dependence* —
  /// i.e., the address being accessed may itself depend on a previously loaded
  /// index, characteristic of indirect memory accesses such as A[B[i]].
  ///
  /// \param Slice - The set of instructions forming the backward slice of a
  /// load address.
  /// \returns true if the slice contains an integer-producing load, false
  /// otherwise.
  bool sliceHasIndexDependence(
      const SmallPtrSet<const Instruction *, 32> &Slice) const {

    // Iterate through all instructions in the slice.
    for (const Instruction *I : Slice) {

      // Focus only on load instructions — potential data dependencies.
      if (const LoadInst *LI = dyn_cast<LoadInst>(I)) {

        // Check if the load produces an integer value (e.g., i = B[j]).
        // This indicates that the current load’s address computation may rely
        // on this value.
        if (LI->getType()->isIntegerTy() &&
            !LI->getPointerOperand()->getType()->isPointerTy()) {
          DBG("    [IMA-check] Slice contains integer load: " << *LI << "\n");

          // Presence of such a load implies an index-based memory access
          // pattern.
          return true;
        }
      }
    }

    // If no integer-producing loads were found, this slice is not
    // index-dependent.
    return false;
  }

  /// \brief Classifies a given LoadInst as Irregular Memory Access (IMA),
  /// PointerChase, Direct, or Other, based on static analysis of its address.
  ///
  /// This routine integrates multiple analysis signals:
  ///  - Whether the load’s pointer originates from a call (e.g., STL vector
  ///  access).
  ///  - Whether it exhibits index-based dependence (via slicing).
  ///  - Whether its access pattern is affine (via ScalarEvolution).
  ///  - Whether it matches known structural forms (e.g., A[B[i]] or pointer
  ///  chasing).
  ///
  /// \param L     The load instruction being analyzed.
  /// \param Slice The backward slice of the load’s pointer operand.
  /// \param LI    LoopInfo, used for loop-based classification.
  /// \param SE    ScalarEvolution, used to check for affine index evolution.
  /// \returns One of {"IMA", "PointerChase", "Direct", "Other"}.
  StringRef
  classifyIrregularAccess(const LoadInst *L,
                          const SmallPtrSet<const Instruction *, 32> &Slice,
                          LoopInfo &LI, ScalarEvolution &SE) {

    const Value *Ptr = L->getPointerOperand();
    DBG("  [Classify] Analyzing: " << *L << "\n");

    // === Case 1: Pointer operand originates from a function call ===
    // This usually represents STL-style data structure access, e.g.:
    //    %idx = call @std::vector::operator[](this, i)
    if (const CallInst *CI = dyn_cast<CallInst>(Ptr)) {
      if (isVectorSubscriptCall(CI)) {
        DBG("    -> Detected STL vector operator[] call: "
            << CI->getCalledFunction()->getName() << "\n");

        // If the call has at least 2 arguments, the second is likely the index.
        if (CI->arg_size() >= 2) {
          const Value *Index = CI->getArgOperand(1);

          // --- Case 1.1: Index is directly loaded from memory (A[B[i]] form)
          // ---
          if (isa<LoadInst>(Index)) {
            DBG("    -> Index argument is a load (integer index): IMA.\n");
            return "IMA";
          }

          // --- Case 1.2: Index expression has integer load dependence ---
          // Recompute slice of the index argument itself.
          SmallPtrSet<const Instruction *, 32> IndexSlice =
              computeStaticSlice(const_cast<Value *>(Index), true);
          if (sliceHasIndexDependence(IndexSlice)) {
            DBG("    -> Slice of index shows integer load dependence: IMA.\n");
            return "IMA";
          }

          // --- Case 1.3: ScalarEvolution says index evolves affinely ---
          // If SCEV can prove index follows affine loop pattern, mark as
          // regular.
          if (const SCEV *S = SE.getSCEV(const_cast<Value *>(Index))) {
            if (const auto *AR = dyn_cast<SCEVAddRecExpr>(S)) {
              const Loop *L = AR->getLoop();
              if (SE.isSCEVable(Index->getType()) && !isa<SCEVUnknown>(S) &&
                  L && SE.hasComputableLoopEvolution(S, L)) {
                DBG("    -> Index proven affine (true IV): Regular.\n");
                return "Other";
              }
            }
          }

          // --- Case 1.4: Fallback ---
          // If no affine pattern is found, assume irregular (non-affine) index.
          DBG("    -> Vector subscript depends on non-affine index: IMA.\n");
          return "IMA";
        }

        // If vector subscript call has no index argument, treat as regular.
        DBG("    -> Vector subscript without index argument.\n");
        return "Other";
      }
    }

    // === Case 2: Standard GEP-style array access ===
    // Handles patterns like load A[i] or load A[B[i]].
    const Value *Base = Ptr->stripPointerCasts();
    const GetElementPtrInst *G = dyn_cast<GetElementPtrInst>(Base);
    if (!G) {
      DBG("    -> No GEP: Direct access.\n");
      return "Direct";
    }

    // === Case 3: Pointer-chasing pattern ===
    // Loads that themselves produce pointers (e.g., linked list dereferences).
    if (L->getType()->isPointerTy()) {
      DBG("    -> Load result is a pointer: PointerChase.\n");
      return "PointerChase";
    }

    // === Case 4: Indirect array access (A[B[i]] form) ===
    // Structural pattern matching using the IMA detection helper.
    if (isIndirectArrayAccess(L)) {
      DBG("    -> Structural pattern detected: IMA.\n");
      return "IMA";
    }

    // === Case 5: Slice-based dependence analysis ===
    // If backward slice contains integer load, index-dependent pattern → IMA.
    if (sliceHasIndexDependence(Slice)) {
      DBG("    -> Slice-based pattern detected: IMA.\n");
      return "IMA";
    }

    // === Default case ===
    // If no irregularity is detected, mark as Other (regular access).
    DBG("    -> Default: Other.\n");
    return "Other";
  }
};

#undef DBG

/// \struct LoopBounds
/// \brief Represents the key structural components of a canonical loop.
///
/// This struct captures essential information extracted from a loop header:
///   - The **induction variable (IV)**: typically an integer PHI node.
///   - The **lower bound**: the incoming value from the preheader.
///   - The **upper bound**: the loop exit condition operand.
///   - The **comparison instruction**: ICmpInst controlling the loop exit.
///
/// Used by analyses such as Scalar Evolution or prefetch heuristics to
/// reason about loop iteration ranges.
struct LoopBounds {
  Value *InductionVar = nullptr; ///< The loop’s induction variable (PHI node)
  Value *LowerBound = nullptr; ///< Initial value of IV (from outside the loop)
  Value *UpperBound = nullptr; ///< Loop termination bound (limit in compare)
  ICmpInst *Compare = nullptr; ///< Comparison instruction in the loop header
};

/// \brief Extracts basic loop bound information (IV, lower, upper) from a
/// loop’s header.
///
/// This function attempts to identify a *canonical loop form* such as:
/// ```
/// for (i = 0; i < N; ++i) { ... }
/// ```
/// by analyzing PHI nodes and conditional branches in the loop header.
/// It returns an optional `LoopBounds` structure if successful.
///
/// \param L The loop to analyze.
/// \returns std::optional<LoopBounds> if identifiable; std::nullopt otherwise.
static std::optional<LoopBounds> getLoopBoundsFromHeader(Loop *L) {
  // Retrieve the loop header basic block.
  BasicBlock *Header = L->getHeader();
  if (!Header)
    return std::nullopt;

  LoopBounds B;

  // === 1. Identify the induction PHI node ===
  //
  // Search through all instructions in the loop header for a PHI node
  // that serves as the loop induction variable (commonly the iterator).
  // Prefer integer-typed PHIs since most loops iterate over integers.
  PHINode *IndPhi = nullptr;
  for (Instruction &I : *Header) {
    if ((IndPhi = dyn_cast<PHINode>(&I))) {
      // Prefer integer induction variables for clarity and analyzability.
      if (IndPhi->getType()->isIntegerTy())
        break;
    }
  }
  if (!IndPhi)
    return std::nullopt; // No suitable induction variable found
  B.InductionVar = IndPhi;

  // === 2. Locate the loop's comparison instruction ===
  //
  // Canonical loops typically have a header terminator like:
  //   br i1 %cmp, label %body, label %exit
  // where %cmp is an ICmpInst comparing the IV and an upper bound.
  auto *BI = dyn_cast<BranchInst>(Header->getTerminator());
  if (!BI || !BI->isConditional())
    return std::nullopt; // No conditional branch → not a loop guard

  // Extract the compare instruction used in the branch condition.
  auto *Cmp = dyn_cast<ICmpInst>(BI->getCondition());
  if (!Cmp)
    return std::nullopt;
  B.Compare = Cmp;

  // === 3. Distinguish induction variable vs. upper bound ===
  //
  // The ICmp operands define the loop’s iteration limit.
  // Example:   %cmp = icmp slt i32 %i, %N
  // Here, %i is the IV, and %N is the upper bound.
  Value *Op0 = Cmp->getOperand(0);
  Value *Op1 = Cmp->getOperand(1);
  if (Op0 == IndPhi)
    B.UpperBound = Op1;
  else if (Op1 == IndPhi)
    B.UpperBound = Op0;
  else
    return std::nullopt; // can't match IV < bound pattern

  // === 4. Identify the lower bound ===
  //
  // The PHI node typically has one incoming edge from outside the loop
  // (initial value) and one from inside (incremented IV).
  // The incoming value from a non-loop block gives us the lower bound.
  for (unsigned I = 0; I < IndPhi->getNumIncomingValues(); ++I) {
    BasicBlock *IncomingBB = IndPhi->getIncomingBlock(I);
    if (!L->contains(IncomingBB)) {
      B.LowerBound = IndPhi->getIncomingValue(I);
      break;
    }
  }

  // Successfully extracted canonical loop bounds.
  return B;
}

/// \brief Attempts to find a common index-producing value among multiple
///        Irregular Memory Access (IMA) load instructions.
///
/// This helper scans through all detected IMA loads and identifies whether
/// they share a *common index source* — for example, the same array or vector
/// used to index multiple dependent loads (e.g., `A[B[i]]`, `C[B[i]]`).
///
/// It supports three major structural patterns:
///   1. **CSR-style (GetElementPtr-based)**: canonical GEPs like `A[idx]`
///   2. **STL vector subscript (CallInst)**: e.g., `std::vector::operator[]()`
///   3. **Fallback heuristic**: when neither pattern applies, use direct loads
///
/// \param IMALoads  The list of LoadInsts classified as IMA.
/// \returns Pointer to the shared Value (index producer), or nullptr if not
/// found.
static Value *
findCommonIndexProducer(const SmallVectorImpl<LoadInst *> &IMALoads) {
  // If there are no IMA loads, nothing to analyze.
  if (IMALoads.empty())
    return nullptr;

  // Will store all extracted index operands from the candidate loads.
  SmallVector<Value *, 8> AllIndices;

  // ---------------------------------------------------------------------------
  // Iterate through all IMA loads and extract their index-producing values.
  // ---------------------------------------------------------------------------
  for (LoadInst *L : IMALoads) {
    if (!L)
      continue;

    // Retrieve the pointer operand used by the load.
    Value *Ptr = L->getPointerOperand()->stripPointerCasts();

    // === Case 1: CSR-style (GEP-based) pattern ===
    //
    // Common in sparse matrix formats like CSR, where we see patterns like:
    //    %idx = load i32, i32* %JA
    //    %val = load double, double* getelementptr(%A, %idx)
    //
    // The GetElementPtr (GEP) instruction encodes the index arithmetic.
    if (auto *GEP = dyn_cast<GetElementPtrInst>(Ptr)) {
      Value *Index = nullptr;
      // The last integer operand of GEP is typically the index
      for (unsigned I = 1; I < GEP->getNumOperands(); ++I) {
        if (GEP->getOperand(I)->getType()->isIntegerTy())
          Index = GEP->getOperand(I);
      }
      if (!Index)
        continue;

      // Strip trivial casts
      Index = Index->stripPointerCasts();
      if (auto *SE = dyn_cast<SExtInst>(Index))
        Index = SE->getOperand(0);
      else if (auto *ZE = dyn_cast<ZExtInst>(Index))
        Index = ZE->getOperand(0);
      else if (auto *TR = dyn_cast<TruncInst>(Index))
        Index = TR->getOperand(0);

      // If index is a load (e.g. JA[j]), that’s our index producer
      if (auto *IdxLoad = dyn_cast<LoadInst>(Index))
        AllIndices.push_back(IdxLoad);
      else
        AllIndices.push_back(Index);
    }

    // === Case 2: STL vector-based pattern (CallInst) ===
    //
    // Detects calls to C++ STL vector subscript operator, e.g.:
    //    %val = call double* @_ZNSt6vectorIdEixEv(%vec, %idx)
    // The second argument is the index operand.
    else if (auto *Call = dyn_cast<CallInst>(Ptr)) {
      if (isVectorSubscriptCall(Call)) {
        // std::vector::operator[](base, index)
        // Extract the index argument from operator[] call.
        Value *Index = Call->getArgOperand(1)->stripPointerCasts();

        // Strip trivial casts (same as above)
        if (auto *SE = dyn_cast<SExtInst>(Index))
          Index = SE->getOperand(0);
        else if (auto *ZE = dyn_cast<ZExtInst>(Index))
          Index = ZE->getOperand(0);
        else if (auto *TR = dyn_cast<TruncInst>(Index))
          Index = TR->getOperand(0);

        // Record the index or its producing load.
        if (auto *IdxLoad = dyn_cast<LoadInst>(Index))
          AllIndices.push_back(IdxLoad);
        else
          AllIndices.push_back(Index);
      }
    }

    // === Case 3: Fallback heuristic ===
    //
    // When no explicit GEP or vector call pattern is found, treat the pointer
    // operand itself as a potential index-producing load.
    else if (auto *LI = dyn_cast<LoadInst>(Ptr)) {
      AllIndices.push_back(LI);
    }
  }

  // ---------------------------------------------------------------------------
  // Post-processing: identify if there’s a single shared index producer.
  // ---------------------------------------------------------------------------

  // If no index candidates were extracted, no common index exists.
  if (AllIndices.empty())
    return nullptr;

  // If there’s only one, that’s trivially the common one
  if (AllIndices.size() == 1)
    return AllIndices.front();

  // Otherwise, check if *all* extracted indices refer to the same Value.
  Value *Common = nullptr;
  for (Value *Idx : AllIndices) {
    if (!Common)
      Common = Idx;
    else if (Common != Idx)
      return nullptr; // mismatch — not a shared index producer
  }

  // All identical — success!
  return Common;
}

/// \brief Inserts prefetch instructions into innermost loops that exhibit
///        irregular memory access (IMA) patterns identified by SlicingAnalysis.
///
/// The pass identifies irregular access chains (A[B[i]]) and, when possible,
/// generates speculative prefetches for future iterations, based on loop bounds
/// and prefetch distance thresholds.
///
/// Key algorithmic steps:
///   1. Identify innermost loops in the function.
///   2. Collect IMA loads (tagged with "prefetch_ima" metadata).
///   3. Derive loop bounds (IV, lower, upper) using
///   `getLoopBoundsFromHeader()`.
///   4. Identify a common index producer among IMAs
///   (`findCommonIndexProducer()`).
///   5. Insert a runtime-guarded prefetch block:
///        - Range and degree checks
///        - Prefetch target address computation
///        - Calls to llvm.prefetch intrinsic
///
/// \param F     Function being processed.
/// \param SA    SlicingAnalysis providing load classification metadata.
/// \param LInfo LoopInfo for iterating over loop hierarchy.
/// \param SE    ScalarEvolution to infer loop bounds and index evolution.
/// \param DL    DataLayout for size/type information.
void insertPrefetchesInLoops(Function &F, SlicingAnalysis &SA, LoopInfo &LInfo,
                             ScalarEvolution &SE, const DataLayout &DL) {
  Module *M = F.getParent();
  LLVMContext &Ctx = M->getContext();

  Type *I64Ty = Type::getInt64Ty(Ctx);

  // ---------------------------------------------------------------------------
  // Helper: Extracts a canonical “base address” from a pointer operand.
  //
  // This function recursively peels through pointer casts, GEPs, global vars,
  // and constant expressions to find a meaningful address root suitable for
  // computing the prefetch target.
  // ---------------------------------------------------------------------------
  auto GetAddressForPointerValue = [&](Value *V) -> Value * {
    if (!V)
      return nullptr;
    Value *S = V->stripPointerCasts();

    // Case 1: If it’s a GEP, follow its base pointer.
    if (auto *G = dyn_cast<GetElementPtrInst>(S)) {
      Value *Base = G->getPointerOperand()->stripPointerCasts();
      if (auto *BaseLoad = dyn_cast<LoadInst>(Base))
        return BaseLoad->getPointerOperand()->stripPointerCasts();
      if (auto *GV = dyn_cast<GlobalVariable>(Base))
        return GV;
      return Base;
    }

    // Case 2: If it’s a load, use the load’s pointer operand.
    if (auto *LI = dyn_cast<LoadInst>(S))
      return LI->getPointerOperand()->stripPointerCasts();

    // Case 3: Handle constant expression GEP (from static arrays).
    if (auto *CE = dyn_cast<ConstantExpr>(S)) {
      if (CE->getOpcode() == Instruction::GetElementPtr) {
        Value *Base = CE->getOperand(0)->stripPointerCasts();
        if (auto *BaseLoad = dyn_cast<LoadInst>(Base))
          return BaseLoad->getPointerOperand()->stripPointerCasts();
        if (auto *GV = dyn_cast<GlobalVariable>(Base))
          return GV;
        return Base;
      }
    }

    // Case 4: Global variable directly used.
    if (auto *GV = dyn_cast<GlobalVariable>(S))
      return GV;

    // Default: return as-is.
    return S;
  };

  // ---------------------------------------------------------------------------
  // Process loops in preorder (outermost → innermost) but handle innermost
  // only.
  // ---------------------------------------------------------------------------
  for (Loop *Top : LInfo.getLoopsInPreorder()) {
    SmallVector<Loop *, 8> Work;
    Work.push_back(Top);
    while (!Work.empty()) {
      Loop *L = Work.pop_back_val();
      for (Loop *Sub : L->getSubLoops())
        Work.push_back(Sub);
      if (!L->getSubLoops().empty())
        continue; // skip non-innermost loops

      // --- Step 1: Collect IMA loads within this loop ---
      SmallVector<LoadInst *, 8> ImaLoads;
      for (BasicBlock *BB : L->blocks())
        for (Instruction &I : *BB)
          if (auto *LI = dyn_cast<LoadInst>(&I))
            if (LI->getMetadata("prefetch_ima"))
              ImaLoads.push_back(LI);

      if (ImaLoads.empty())
        continue; // nothing to prefetch here

      // --- Cost Model: Skip if benefit is below threshold ---
      if (EnableCostModel) {
        int benefit = computePrefetchBenefit(L, ImaLoads.size(), SE);
        if (benefit < CostThreshold) {
          // Skip this loop - prefetching not worth the overhead
          continue;
        }
      }

      // --- Step 2: Extract canonical loop bounds ---
      auto BoundsOpt = getLoopBoundsFromHeader(L);
      if (!BoundsOpt)
        continue;
      LoopBounds Bnds = *BoundsOpt;
      Value *LowerBound = Bnds.LowerBound;
      Value *UpperBound = Bnds.UpperBound;
      Value *IndVar = Bnds.InductionVar;
      if (!LowerBound || !UpperBound || !IndVar)
        continue;

      // --- Step 3: Find common index producer across IMA loads ---
      Value *CommonIdxProd = findCommonIndexProducer(ImaLoads);
      if (!CommonIdxProd)
        continue;
      LoadInst *CommonLoad = dyn_cast<LoadInst>(CommonIdxProd);
      if (!CommonLoad)
        continue;

      // Determine insertion point: right after the common index load.
      Instruction *InsertPt = CommonLoad->getNextNode();
      if (!InsertPt)
        InsertPt = CommonLoad->getParent()->getTerminator();

      // --- Step 4: Construct prefetch condition and setup ---
      IRBuilder<> Bldr(InsertPt);
      Value *IV64 = Bldr.CreateSExtOrTrunc(IndVar, I64Ty, "iv64");
      Value *DegThreshVal = ConstantInt::get(I64Ty, (int64_t)DegreeThreshold);

      int effectiveDistance = PrefetchDistance;
      if (AdaptiveDistance) {
        effectiveDistance = computeAdaptiveDistance(L, MemoryLatencyCycles);
      }
      Value *PrefDistVal = ConstantInt::get(I64Ty, (int64_t)effectiveDistance);

      // Degree condition: (upper - lower) > DegreeThreshold
      Value *Diff = Bldr.CreateSub(UpperBound, LowerBound, "deg.diff");
      Value *CmpDeg = Bldr.CreateICmpSGT(Diff, DegThreshVal, "deg.cond");

      // Range condition: (iv + D) < upper
      Value *IVPlusD = Bldr.CreateAdd(IV64, PrefDistVal, "iv.plus.d");
      Value *CmpRange = Bldr.CreateICmpSLT(IVPlusD, UpperBound, "range.cond");

      // Combined predicate for prefetch eligibility
      Value *PrefetchCond = Bldr.CreateAnd(CmpDeg, CmpRange, "prefetch.cond");

      // Create a “then” block guarded by PrefetchCond
      Instruction *ThenTerm =
          SplitBlockAndInsertIfThen(PrefetchCond, InsertPt, false);
      BasicBlock *ThenBlock = ThenTerm->getParent();
      ThenBlock->setName("prefetch.then");
      BasicBlock *ContBlock = ThenBlock->getSingleSuccessor();
      if (ContBlock)
        ContBlock->setName("prefetch.cont");
      IRBuilder<> ThenB(ThenTerm);

      // -----------------------------------------------------------------------
      // Step 5: Derive JA[iv + d] — the "next source" index for prefetch target
      // -----------------------------------------------------------------------
      Value *JaLoadPtr = CommonLoad->getPointerOperand();
      Value *JaLoadPtrStripped = JaLoadPtr->stripPointerCasts();
      Type *JAElemTy = CommonLoad->getType(); // usually i32

      Value *NextSrcI64 = nullptr;

      // Case 1: std::vector::operator[]
      if (auto *Call = dyn_cast<CallInst>(JaLoadPtrStripped)) {
        if (isVectorSubscriptCall(Call)) {
          Value *FutureIdx = ThenB.CreateSExtOrTrunc(
              IVPlusD, Call->getArgOperand(1)->getType(), "future.idx");
          // Perform call: JA[iv + D]
          Value *FutureElemPtr = ThenB.CreateCall(
              Call->getCalledFunction(), {Call->getArgOperand(0), FutureIdx},
              "ja.future.ptr");
          Value *FutureVal =
              ThenB.CreateLoad(JAElemTy, FutureElemPtr, "ja.future.val");
          NextSrcI64 =
              ThenB.CreateSExtOrTrunc(FutureVal, I64Ty, "next_src.i64");
        }
      }

      // Case 2: Normal GEP-based
      if (!NextSrcI64) {
        Value *JaReloadAddr = GetAddressForPointerValue(JaLoadPtrStripped);
        if (!JaReloadAddr || !JaReloadAddr->getType()->isPointerTy())
          continue;
        Type *JaLoadValueType = PointerType::get(JAElemTy, 0);
        Value *JaPtrLocal =
            ThenB.CreateLoad(JaLoadValueType, JaReloadAddr, "ja.ptr.loaded");
        Value *IdxForJA = ThenB.CreateSExtOrTrunc(IVPlusD, I64Ty, "idx_for_ja");
        Value *JaGep =
            ThenB.CreateGEP(JAElemTy, JaPtrLocal, IdxForJA, "ja.next.gep");
        Value *NextSrc = ThenB.CreateLoad(JAElemTy, JaGep, "next_src");
        NextSrcI64 = ThenB.CreateSExtOrTrunc(NextSrc, I64Ty, "next_src.i64");
      }

      if (!NextSrcI64)
        continue;

      // -----------------------------------------------------------------------
      // Step 6: Insert llvm.prefetch() calls for dependent irregular arrays
      // -----------------------------------------------------------------------
      for (LoadInst *LI : ImaLoads) {
        if (LI == CommonLoad)
          continue;

        Value *PtrOp = LI->getPointerOperand()->stripPointerCasts();
        Type *ElemTy = LI->getType();

        // Case A: std::vector::operator[]
        if (auto *Call = dyn_cast<CallInst>(PtrOp)) {
          if (isVectorSubscriptCall(Call)) {
            // Create future pointer = Vec[idx + PrefetchDistance]
            Value *FutureIdx = ThenB.CreateSExtOrTrunc(
                NextSrcI64, Call->getArgOperand(1)->getType(), "future.idx");

            Value *FuturePtr = ThenB.CreateCall(
                Call->getCalledFunction(), {Call->getArgOperand(0), FutureIdx},
                "future.ptr");

            // Bitcast to i8* (as required by llvm.prefetch intrinsic)
            Type *I8PtrTy = PointerType::get(ThenB.getInt8Ty(), 0);
            Value *FuturePtrI8 =
                ThenB.CreateBitCast(FuturePtr, I8PtrTy, "future.i8");

            // Prefetch parameters: RW=0(read), locality=LocalityHint,
            // cachetype=1(data)
            Value *RW = ThenB.getInt32(0);
            Value *Locality = ThenB.getInt32(LocalityHint);
            Value *CacheType = ThenB.getInt32(1);

            ThenB.CreateIntrinsic(Intrinsic::prefetch, {FuturePtrI8->getType()},
                                  {FuturePtrI8, RW, Locality, CacheType});

            // errs() << "[PrefetchIMA] Inserted vector prefetch for: " << *LI
            //        << "\n";
            continue;
          }
        }

        // Case B: Normal GEP/global pointer access - with batch prefetching
        Value *BaseAddr = GetAddressForPointerValue(PtrOp);
        if (!BaseAddr || !BaseAddr->getType()->isPointerTy())
          continue;

        Type *LoadValueTy = PointerType::get(ElemTy, 0);
        Value *BasePtr =
            ThenB.CreateLoad(LoadValueTy, BaseAddr, "base.ptr.loaded");

        // Batch prefetching: prefetch BatchSize consecutive addresses
        // This leverages spatial locality from graph reordering
        int effectiveBatchSize = std::min(std::max(BatchSize.getValue(), 1), 8);
        for (int batch = 0; batch < effectiveBatchSize; batch++) {
          Value *BatchOffset = ConstantInt::get(I64Ty, batch);
          Value *BatchIdx =
              ThenB.CreateAdd(NextSrcI64, BatchOffset, "batch.idx");

          Value *Gep =
              ThenB.CreateGEP(ElemTy, BasePtr, BatchIdx, "prefetch.gep");

          // Cast to i8* for llvm.prefetch
          Type *I8PtrTy = PointerType::get(ThenB.getInt8Ty(), 0);
          Value *GepI8 = ThenB.CreateBitCast(Gep, I8PtrTy, "prefetch.i8");

          Value *RW = ThenB.getInt32(0);
          Value *Locality = ThenB.getInt32(LocalityHint);
          Value *CacheType = ThenB.getInt32(1);

          ThenB.CreateIntrinsic(Intrinsic::prefetch, {GepI8->getType()},
                                {GepI8, RW, Locality, CacheType});
        }
      }
    }
  }
}

} // namespace

class PrefetchIMAPass : public PassInfoMixin<PrefetchIMAPass> {
public:
  PreservedAnalyses run(Function &F, FunctionAnalysisManager &AM) {
    // Check function pattern: if specified, match against it; otherwise process
    // all
    if (!FunctionPattern.empty()) {
      std::string FnName = F.getName().str();
      // Simple substring match
      if (FnName.find(FunctionPattern) == std::string::npos) {
        return PreservedAnalyses::all();
      }
    }
    // Note: If FunctionPattern is empty, we process all functions that have IMA
    // loads

    auto &LInfo = AM.getResult<LoopAnalysis>(F);
    auto &SE = AM.getResult<ScalarEvolutionAnalysis>(F);

    SlicingAnalysis Slicing(F, AM, LInfo, SE);

    insertPrefetchesInLoops(F, Slicing, LInfo, SE,
                            F.getParent()->getDataLayout());

    return PreservedAnalyses::all();
  };
};

extern "C" LLVM_ATTRIBUTE_WEAK ::llvm::PassPluginLibraryInfo
llvmGetPassPluginInfo() {
  return {.APIVersion = LLVM_PLUGIN_API_VERSION,
          .PluginName = "Prefetching Irregular Memory Access Pass",
          .PluginVersion = "v0.1",
          .RegisterPassBuilderCallbacks = [](PassBuilder &PB) {
            PB.registerPipelineParsingCallback(
                [](StringRef Name, FunctionPassManager &FPM,
                   ArrayRef<PassBuilder::PipelineElement>) {
                  if (Name == "prefetch-ima") {
                    FPM.addPass(PrefetchIMAPass());
                    return true;
                  }
                  return false;
                });
          }};
}
