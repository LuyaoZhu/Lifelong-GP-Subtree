# LEGO: Lifelong Learning in Genetic Programming with Building Block Consolidation for Dynamic Flexible Job Shop Scheduling

This repository provides the implementation of **LEGO** (Lifelong GP with building block consolidation), a subtree-level lifelong learning framework for the automated design of scheduling rules in Dynamic Flexible Job Shop Scheduling (DFJSS), as proposed in the submitted manuscript.

This repository is developed based on an open-source Genetic Programming for Job Shop Scheduling (GPJSS) benchmarking framework built upon the Java ECJ library. The original package structure and foundational scheduling simulation modules have been retained. On top of the multi-tree representation, continuous task transitions and surrogate-assisted lifelong learning mechanisms (LLSGP), this repository introduces the subtree-level building block consolidation mechanisms proposed in the manuscript.

## ✨ Overview

Existing lifelong GP methods mitigate catastrophic forgetting only at the individual level. As a result, standard genetic operators can easily destroy the task-critical subtrees (building blocks) inside GP individuals. LEGO addresses this with two strategies:

- **Adaptive Parent Selection.** GP individuals are grouped into three knowledge groups: *historical generalists*, *current specialists* and *current generalists*. A multi-armed bandit (UCB with Softmax selection) adaptively determines how often each crossover parent combination is used.
- **Building Block Consolidation.** Subtree importance is measured by how much a parent's decision-making behaviour changes when the subtree is replaced by its mean output. Behaviour change is quantified by the Spearman rank correlation of candidate rankings over sampled decision points. These task-aware importance values then guide crossover, so that important subtrees are preferentially preserved and donated, while unimportant ones are replaced.

## 🛠️ Environment Setup & Dependencies

**Prerequisites:**
- Java JDK 11 or higher.
- IntelliJ IDEA (recommended project environment).

**Project Structure Realignment:**
- Make sure that the core algorithmic source files in the `src/` folder and the external library dependencies (e.g., Apache Commons Math, ECJ core files) are fully loaded in your IDE project build path.
- If you use IntelliJ IDEA, set the project root directory as the **Content Root** and mark the `src/` folder as **Sources Root**.

**Code Modules:**
- `src/ec/`: Core ECJ foundational framework functions.
- `src/yimei/jss/`: Domain-specific packages for Job Shop Scheduling, containing dispatching rule representations, discrete event simulation environments and tree evaluation models.

## 📂 Key Architecture Modifications

The core components of LEGO are located in `src/yimei/jss/algorithm/lifelongGP-subtree/`:

- **Individual group management:** collection and maintenance of historical generalists (HG), current specialists (CS) and current generalists (CG) across the two learning stages of each task.
- **Multi-armed bandit parent selection:** UCB-based reward estimation for the four crossover parent combinations (HG+CS, HG+CG, CS+CG, CG+CG), with Softmax-based probabilistic arm selection.
- **Subtree importance evaluation:** replacement-based behaviour change measurement on sampled decision points of each seen task, with group-specific importance assignment.
- **Importance-guided crossover:** a crossover operator that selects the removed and donated subtrees according to their importance values.

Mechanisms inherited from the surrogate-assisted lifelong GP framework (two-stage learning, surrogate model construction and generalist preselection) are reused without modification.

## 🚀 Running Experiments: Launching LEGO

1. Locate the target parameter configuration file:
   ```
   src/yimei/jss/algorithm/lifelongGP-subtree/multipletreegp-dynamic.params
   ```
2. Run the main class `src/yimei/jss/gp/GPRun` with the following command-line argument:
   ```
   -file src/yimei/jss/algorithm/lifelongGP-subtree/multipletreegp-dynamic.params
   ```

## ⚙️ Default Parameter Settings

| Parameter | Value |
|---|---|
| Number of tasks | 4 |
| Generations per task | 100 |
| Population size | 500 |
| Crossover / Mutation / Reproduction rate | 80% / 15% / 5% |
| HG / CS archive ratio | 5% |
| UCB exploration coefficient *c* | 1 |
| Softmax temperature τ | 0.2 |
| Decision points per task for importance evaluation | 100 |

Other GP settings (tournament size, tree depth limits, initialisation method, etc.) follow the parameter file above and the manuscript.
