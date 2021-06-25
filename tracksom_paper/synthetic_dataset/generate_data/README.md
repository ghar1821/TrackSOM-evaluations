# synthetic dataset generator

This package contains all modules used to generate synthetic dataset using gaussian distribution.

The actual dataset is generated for the ChronoClust paper
```
@article{putri2019chronoclust,
  title={ChronoClust: Density-based clustering and cluster tracking in high-dimensional time-series data},
  author={Putri, Givanna H and Read, Mark N and Koprinska, Irena and Singh, Deeksha and R{\"o}hm, Uwe and Ashhurst, Thomas M and King, Nicholas JC},
  journal={Knowledge-Based Systems},
  volume={174},
  pages={9--26},
  year={2019},
  publisher={Elsevier}
}
```

The file `dataset_generator.py` and `distributions_generators.py` are used by the ChronoClust paper's experiment to generate the synthetic dataset.

The code `extract_data_without_noise.R` is used to remove the 100 noise data points.

The code `add_disappearing_cluster.py` is used to add a suddenly appearing and disappearing cluster in day 3-5.
