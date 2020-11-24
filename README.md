# Inverse LQR with automatic differentiation

This is a simple example that accompanies the technical note [Automatic differentiation of Sylvester, Lyapunov, and algebraic Riccati equations](https://arxiv.org/abs/2011.11430).

# Run example
```sh
mkdir results
dune exec ./inverse_lqr.exe -- -d results
```
# Dependencies
- [owl](https://github.com/owlbarn/owl.git)
- [base](https://github.com/janestreet/base)
- [stdio](https://github.com/janestreet/stdio)
- [cmdargs](https://github.com/hennequin-lab/cmdargs)

