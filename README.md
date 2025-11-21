# certifiable-rwhe-calibration
Certifiably globally optimal generalized robot-world and hand-eye calibration. Our algorithm is the first extrinsic calibration method that can jointly solve for the poses of multiple sensors *and* targets, consider an unknown target scale, and provide a computational certificate of global optimality for its maximum likelihood objective function. 

## Usage
See `experiments/rw_multi_eye_multi_hand.jl` for an example of how to use our software. *A more convenient interface is coming soon!*


##  Citation
If you use this work in your research, please cite the [following paper](https://arxiv.org/abs/2507.23045):

```bibtex
@article{wise2025certifably,
  title={A Certifably Correct Algorithm for Generalized Robot-World and Hand-Eye Calibration},
  author={Wise, Emmett and Kaveti, Pushyami and Chen, Qilong and Wang, Wenhao and Singh, Hanumant and Kelly, Jonathan and Rosen, David M and Giamou, Matthew},
  journal={arXiv preprint arXiv:2507.23045},
  year={2025}
}
```