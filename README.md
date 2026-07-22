# Moderna Challenge for WISER Summer Program 2026
author: Sergey Grigorovich

### 1. Introduction and problem statement

mRNA can fold on itself when complementary nucleotides form base pairs, with the resulting secondary structure influencing stability, translation, and manufacturing. Predicting these structures is relevant to mRNA study and design.

Classical tools such as ViennaRNA predict low-energy RNA structures using thermodynamic models. In this project, ViennaRNA provides a reference structure and energy for quantum model comparison. The goal is not to replace this approach, but to study how RNA folding can be represented in a quantum optimization problem.

#### IBM-Moderna studies of mRNA folding (selected)
Recent IBM–Moderna research provides a clear progression toward larger and more hybrid quantum workflows.
- *Alevras et al. (2024)* represented RNA folding as a binary optimization problem and executed a variational algorithm on IBM quantum processors. It is a clear demonstation that larger RNA optimization models could be studied on current hardware rather than only through small simulations.

- *Kumar, Alevras, Metkar et al. (2025)* extended this direction through a broader hybrid workflow. The approach combined quantum execution with classical transformations, and local optimization to address larger RNA sequences.

- *Friedhoff, Metkar, Davis, Kumar, and Galda (2026)* focused on two remaining barriers: the number of qubits required and the difficulty of decoding dense, highly constrained optimization problems. Their compressed encoding and problem-aware decoder moved responsibility into a classical postprocessing stage, showing that the representation of constraints is as important as the algorithm itself.

#### The present project 
The present project continues the direction of these studies by asking a question about the role of structural constraints. Encoding **every** constraint directly in the quantum objective may produce more valid structures, but it adds multiple interactions and can make circuits **more demanding** by increasing depth. Leaving some constraints outside can simplify the circuit, but then invalid selections must be corrected classically.

The project asks a question about the **trade-off** through three versions of the RNA stem-selection problem:
- the **strict** variant places both nucleotide-overlap and crossing-stem constraints inside the QUBO;
- the **relaxed** variant keeps overlap constraints in the QUBO but handles crossings after sampling;
- the **postprocessed** variant uses the simplest quantum objective and leaves both conflict types to classical repair.

All three versions use the same RNA sequences, candidate stems, stem rewards, quantum workflow, repair logic, and ViennaRNA evaluation. The main difference is **where** structural conflicts are resolved. This design allows the project to compare QUBO connectivity, quantum circuit requirements, raw structural validity, dependence on classical repair, and the quality of structures present in the sampled candidate set.

The quantum objective used in this project is a *structural stem-selection proxy*, rather than a thermodynamic Minimum Free Energy (MFE) model. Each candidate stem receives a reward equal to its length. ViennaRNA energy is not included in the QUBO or used during QAOA parameter optimization. It is applied after sampling to evaluate the repaired structures against a classical reference. 

This project is not expected to purpose a better or universal replacement for classical RNA-folding software. It explores **whether stricter quantum constraint encoding provides enough practical benefit to justify its additional circuit cost, or whether a simpler hybrid strategy offers a better balance for current quantum hardware**.

### 2. Data, simulations and hardware

The project relies on synthetic segments of RNA sequences, produced by concatenating real sequences from the BEACON dataset (*Ren, Yuchen, et al., 2024*) noncoding-RNA task and dividing them into target lengths. This approach preserves some balance between biological reality and the need for short, fixed-lenght sequences for testing and scaling comparisons.

Сompleted 1,170 successful variant–sequence runs: **630 Aer simulator** runs and **540 IBM hardware** runs:
- The simulator experiment covered **30 sequences** at each of 7 lengths `10, 12, 14, 16, 18, 19, 20` nucleotides, with all three variants.
- The hardware experiment covered **15 sequences** at each of 12 lengths `10, 15, 20, 25, 30, 35, 37, 40, 41, 42, 43, 44` nucleotides, with all three variants. 

Simulation and optimization were performed on Google Colab virtual machine (Intel Xeon CPU ~2.20 GHz). Simulation runs took **~14 minutes real-time** execution.

Hardware runs were executed circuits with up to **143** logical qubits on the **156-qubit** `ibm_quebec` backend. Recorded **QPU usage** reached **~28 minutes**. However, with transpilation, data transfer, results retrieval, with no other pendig jobs and queues, total execution time reached **~2.2 hours real-time**.

<img src="figures/01_candidate_stem_growth.png" width="400" height="300"> <img src="figures/08_num_qubits.png" width="400" height="300">

The binary search space is defined by the enumerated candidate stems. Candidate stems must contain at least two consecutive allowed base pairs and satisfy the configured minimum loop length. Final structures are therefore limited to noncrossing combinations of stems present in this candidate set.

The representation may not contain every base pair used by the ViennaRNA MFE structure. It does not directly represent isolated single base pairs, pseudoknots, or structural elements that cannot be assembled from the enumerated consecutive stems. As a result, the ViennaRNA reference may be unreachable even with an exact optimizer.

A *representation oracle* is the best valid structure that can be constructed from the complete candidate-stem set. Its difference from the ViennaRNA reference measures **error introduced by the representation itself**.

A sampled oracle is the best repaired structure among the bitstrings actually observed in one solver run. Its difference from the representation oracle also includes sampling and optimization limitations.


### 3. Results

The experiments show a clear trade-off between quantum constraint encoding and classical repair.
Compared with *strict* encoding, the *relaxed* variant **reduced the mean number of quadratic interactions by 37%** on average, across various sequence lenghts, and **reduced mean circuit depth by 21.3%**. The paired comparison showed lower *relaxed*-circuit depth at every tested nontrivial sequence length. See figures below:

<img src="figures/02a_aer_qubo_interaction_growth.png" width="400" height="300"> <img src="figures/02b_ibm_qubo_interaction_growth.png" width="400" height="300">
<img src="figures/07_circuit_depth.png" width="400" height="300"> <img src="figures/12_circuit_depth_relaxed_minus_strict.png" width="400" height="300">

This reduction **did not worsen** repaired the **quality** of the best repaired candidate substantially. Using the *sampled-oracle selection described above*, the mean ViennaRNA energy gap was **0.443 kcal/mol** for *relaxed* and **0.475 kcal/mol** for *strict*. Mean base-pair F1 was **0.801 and 0.811**, respectively. See figures below:

<img src="figures/12_energy_gap_relaxed_minus_strict.png" width="400" height="300"> <img src="figures/12_pair_f1_relaxed_minus_strict.png" width="400" height="300">

**Raw validity** provides a separate view that does not depend on ViennaRNA-based candidate selection. *Strict* encoding improved raw validity on shorter instances because both overlap and crossing conflicts were included in the QUBO. This advantage decreased with problem size. From length 35 onward, the probability-weighted raw validity of both *strict* and *relaxed* hardware samples was close to zero, and the final structures depended strongly on classical repair. See figures below:

<img src="figures/03_raw_validity.png" width="400" height="300"> <img src="figures/12_probability_weighted_raw_validity_relaxed_minus_strict.png" width="400" height="300">

The *postprocessed* variant **minimized circuit cost**, with depth remaining near 6, but lead to **much more work for classical repair**. Because this QUBO contains no conflict penalties, its unconstrained optimum is to select all candidate stems. The quantum objective itself does not distinguish valid from conflicting stem combinations in this variant.

It required **48.91** stem removals on average, compared with **22.50** for *relaxed*, a 117% increase. Mean energy gap vs *relaxed* was also about **3.5 times larger**, **1.546** versus **0.443** kcal/mol, while mean F1 fell from **0.801** to **0.571**. See figures below:

<img src="figures/04_repair_burden.png" width="400" height="300"> <img src="figures/12_probability_weighted_stems_removed_postprocessed_minus_relaxed.png" width="400" height="300">
<img src="figures/05_energy_gap.png" width="400" height="300"> <img src="figures/12_energy_gap_postprocessed_minus_relaxed.png" width="400" height="300">

For more details, check <a href="rna_qubo_results_analysis.ipynb">**analysis notebook**</a> and individual <a href="figures">**figures**</a>

The current results support a comparison of **constraint placement, sampled candidate availability, repair burden, and circuit resources**.

### 4. Summary

In this project, synthetic RNA sequences were tested on quantum hardware up to a length of 44 nucleotides, with selection bias toward sequences having manageable candidate-stem counts, in three modes, defined by the stage of resolving structural constraints relative to the quantum encoding and computation. 

The strongest result is the **observed resource and constraint-placement trade-off**.

Overall, the results of simulations and hardware runs demonstrate that the *relaxed* variant provided **reasonably good balance** among the tested strategies. It substantially **reduced QUBO connectivity and circuit depth** compared to the *strict* encoding while **preserving similar repaired quality**. The *postprocessed* mode minimizes quantum computation complexity, but puts all the burden on classical postprocessing.

The *strict* and *relaxed* variants had similar **sampled-oracle quality**, meaning that they produced repaired candidates with similar best ViennaRNA energy gaps and base-pair F1 values. Because the reference energy and structure were used to select these candidates after sampling, this comparison **should not** be interpreted as native top-1 prediction accuracy. 

*Separate* highest-probability, lowest-QUBO-objective analysis is required to determine how much of the final quality comes from the quantum distribution, the candidate representation, and the classical decoder.

Based on the metrics trend with increasing RNA length, it is likely but uncertain, whether it can be extrapolated to high qubit count and longer sequences. Also, sequence length alone is not sufficient to describe optimization difficulty, because the number of candidate stems and quadratic interactions depends on sequence composition.


*placeholder: link to presentation/video*

### 5. Project workflow

*placeholder: image of workflow chart*

- The workflow begins by loading BEACON source sequences and generating synthetic ones from segments of fixed lengths.
- ViennaRNA is used to calculate an MFE reference. 
- Candidate stems are sorted, less complex variants are selected, enumerated, conflicts are identified, and strict, relaxed, and postprocessed QUBOs are constructed.
- QAOA parameters are optimized on Aer for the simulation experiment. 
- Run the simulation experiments and collect bitstrings.
- Hardware runs use fixed parameters obtained either through full Aer optimization (for sequences 10-20 length) or transferred from the simulated sequence of the same variant with max number of variables (for 20+ nucleotides).
- All measured bitstrings are decoded, invalid stems are repaired
- Each repaired structure is evaluated with ViennaRNA. The current energy-gap and F1 figures use a sampled-oracle summary that selects the best observed candidate using reference information.
- Each repaired structure is evaluated on validity, repair burden, QUBO complexity, circuit resources and runtime.
- Execution results are saved as tables and processed separately in the analysis notebook.

#### QUBO objective

For each candidate stem $i$:

- $x_i \in \{0,1\}$ indicates whether stem $i$ is selected;
- $l_i$ is the number of base pairs in stem $i$;
- $O$ is the set of nucleotide-overlap conflicts;
- $C$ is the set of crossing-stem conflicts.

The **strict** QUBO objective is:

```math
Q_{\mathrm{strict}}(x)
=
-\sum_i l_i x_i
+
P_{\mathrm{overlap}}
\sum_{(i,j)\in O} x_i x_j
+
P_{\mathrm{crossing}}
\sum_{(i,j)\in C} x_i x_j
```

The **relaxed** variant excludes crossing conflicts from the quantum objective:

```math
Q_{\mathrm{relaxed}}(x)
=
-\sum_i l_i x_i
+
P_{\mathrm{overlap}}
\sum_{(i,j)\in O} x_i x_j
```

The **postprocessed** variant contains only the stem-selection rewards:

```math
Q_{\mathrm{postprocessed}}(x)
=
-\sum_i l_i x_i
```

The implementation uses the following penalties:

```math
P_{\mathrm{overlap}}
=
P_{\mathrm{crossing}}
=
1+\max_i l_i
```

Because the penalty is larger than the reward of any individual stem, removing one stem from an encoded conflicting pair always lowers the objective. Therefore, an exact optimum of the strict QUBO cannot contain an encoded overlap or crossing conflict, while an exact optimum of the relaxed QUBO cannot contain an encoded overlap conflict.

These penalties enforce structural constraints and should not be interpreted as thermodynamic free-energy parameters. Penalty sensitivity was not evaluated in the present study.

For more details, check <a href="rna_qubo_execution.ipynb">**execution notebook**</a>.

#### Modules

The notebooks were executed in a hosted Google Colab runtime 1.0.0, with Python 3.12.13 and CPU-only.
<a href="colab-runtime-info.json">**Runtime report**</a>
<a href="requirements-colab-freeze.txt">**Exact package snapshot**</a>

For installation:
`python -m pip install -r requirements.txt`

`/src` modules provide functions for execution and analysis notebooks:

<a href="src/data.py">`data.py`</a> : loading BEACON data, generation of fixed-length synthetic sequences, ViennaRNA reference structures and energies, preparation of the processed sequence table

<a href="src/model.py">`model.py`</a>: enumeration of candidate stems, overlap and crossing conflicts, construction of the *strict*, *relaxed*, and *postprocessed* QUBO versions

<a href="src/quantum.py">`quantum.py`</a>: solver, Aer QAOA optimization and sampling, IBM backend preparation and sampling

<a href="src/analysis.py">`analysis.py`</a>: decoding solver outputs, structural repair, evaluation of repaired structures, aggregation of results, plotting and summary utilities



### 6. Limitations

- Sequences are synthetic segments created from concatenated BEACON sequences
- MFE reference is based on ViennaRNA, rather than experimental ground truth
- Basic hardware-aware optimization
- No advanced encoding/transpilation
- No quantum error correction techniques
- Noise-less simulations
- The QUBO rewards stem length and does not implement a thermodynamic energy model
- Internal loops and other structural features are not represented directly
- Predictions do not include pseudoknots
- Sampled-oracle quality is an upper bound on the observed candidates and is not equivalent to the quality of the highest-probability QAOA output
- Shallow QAOA with p = 1 and limited classical optimization budgets
- Transferred parameters for longer sequences for hardware runs
- Biased sequence selection for manageable candidate-stem counts
- Hardware results come from one backend and limited repetitions, no generalization

### References
1. Alevras, Dimitris, et al. "mRNA secondary structure prediction using utility-scale quantum computers." 2024 IEEE International Conference on Quantum Computing and Engineering (QCE). Vol. 1. IEEE, 2024.
2. Kumar, Vaibhaw, et al. "Towards secondary structure prediction of longer mrna sequences using a quantum-centric optimization scheme." 2025 IEEE International Conference on Quantum Computing and Engineering (QCE). Vol. 1. IEEE, 2025.
3. Friedhoff, Triet, et al. "Pauli Correlation Encoding for mRNA Secondary Structure Prediction: Problem-Aware Decoding for Dense-Constraint QUBOs." arXiv preprint arXiv:2605.20163 (2026).
4. Ren, Yuchen, et al. "Beacon: Benchmark for comprehensive rna tasks and language models." Advances in Neural Information Processing Systems 37 (2024): 92891-92921.

#### Disclosure - Generative AI usage

*Gemini* 2.5 in Google Colab environment - syntax autocompletion, troubleshooting, debugging
*ChatGPT* 5.5 - debugging, technical writing (comments, markdown text cells), editing (grammar, style), pre-release critique