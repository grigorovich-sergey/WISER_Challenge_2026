# WISER-Moderna Challenge for WISER Summer Program 2026
author

### 1. Introduction and problem statement

mRNA can fold on itself when complementary nucleotides form base pairs, with the resulting secondary structure influencing stability, translation, and manufacturing. Predicting these structures is relevant to mRNA study and design.

Classical tools such as ViennaRNA predict low-energy RNA structures using thermodynamic models. In this project, ViennaRNA provides a reference structure and energy for quantum model comparison. The goal is not to replace this approach, but to study how RNA folding can be represented in a quantum optimization problem.

#### IBM-Moderna studies of mRNA folding (selected)
Recent IBM–Moderna research provides a clear progression toward larger and more hybrid quantum workflows.
- Alevras et al. (2024) represented RNA folding as a binary optimization problem and executed a variational algorithm on IBM quantum processors. It is a clear demonstation that larger RNA optimization models could be studied on current hardware rather than only through small simulations.

- Kumar, Alevras, Metkar et al. (2025) extended this direction through a broader hybrid workflow. The approach combined quantum execution with classical transformations, and local optimization to address larger RNA sequences.

- Friedhoff, Metkar, Davis, Kumar, and Galda (2026, *Pauli Correlation Encoding for mRNA Secondary Structure Prediction: Problem-Aware Decoding for Dense-Constraint QUBOs*) focused on two remaining barriers: the number of qubits required and the difficulty of decoding dense, highly constrained optimization problems. Their compressed encoding and problem-aware decoder moved responsibility into a classical postprocessing stage, showing that the representation of constraints is as important as the algorithm itself.

#### The present project 
The present project continues the direction of these studies by asking a question about the role of structural constraints. **Encoding every constraint** directly in the quantum objective may produce more valid structures, but it adds multiple interactions and can make circuits **more demanding** by increasing depth. Leaving some constraints outside can simplify the circuit, but then invalid selections must be corrected classically.

The project asks a question about the **trade-off** through three versions of the RNA stem-selection problem:
- the **strict** variant places both nucleotide-overlap and crossing-stem constraints inside the QUBO;
- the **relaxed** variant keeps overlap constraints in the QUBO but handles crossings after sampling;
- the **postprocessed** variant uses the simplest quantum objective and leaves both conflict types to classical repair.

All three versions use the same RNA sequences, candidate stems, stem rewards, quantum workflow, repair logic, and ViennaRNA evaluation. The main difference is **where** structural conflicts are resolved. This design allows the project to compare quantum resource requirements, raw structural validity, dependence on classical repair, and agreement with the ViennaRNA.

This project is not expected to purpose a better or universal replacement for classical RNA-folding software. It explores *whether stricter quantum constraint encoding provides enough practical benefit to justify its additional circuit cost, or whether a simpler hybrid strategy offers a better balance for current quantum hardware*.

### 2. Data, simulations and hardware

The project relies on synthetic segments of RNA sequences, produced by concatenating real sequences from the BEACON dataset (Ren, Yuchen, et al., 2024) noncoding-RNA task and dividing them into target lengths. This approach preserves some balance between biological reality and the need for short, fixed-lenght sequences for testing and scaling comparisons.
*placeholder: lenghts and numbers for simulated and hardware runs*

*placeholder: simulator and hardware runs: hardware specifications, real time used, counts, QPU time*
*placeholder: graphs for qubit per lenght, variables per lenght, compute time per lenght*

### 3. Results
*placeholder: brief summary description, key findings*
*placeholder: detailed results + graphs links*
*placeholder: link to analyzis notebook*
*placeholder: link to presentation/video*

### 4. Project workflow
*placeholder: general overview*
*placeholder: link to execution notebook*
*placeholder: image of workflow chart*
*placeholder: short modules description*


### 5. References
1. Alevras, Dimitris, et al. "mRNA secondary structure prediction using utility-scale quantum computers." 2024 IEEE International Conference on Quantum Computing and Engineering (QCE). Vol. 1. IEEE, 2024.
2. Kumar, Vaibhaw, et al. "Towards secondary structure prediction of longer mrna sequences using a quantum-centric optimization scheme." 2025 IEEE International Conference on Quantum Computing and Engineering (QCE). Vol. 1. IEEE, 2025.
3. Friedhoff, Triet, et al. "Pauli Correlation Encoding for mRNA Secondary Structure Prediction: Problem-Aware Decoding for Dense-Constraint QUBOs." arXiv preprint arXiv:2605.20163 (2026).
4. Ren, Yuchen, et al. "Beacon: Benchmark for comprehensive rna tasks and language models." Advances in Neural Information Processing Systems 37 (2024): 92891-92921.