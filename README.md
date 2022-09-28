# eigenvalue-problem

## problem:

$$
\begin{cases}
\frac{y_{n+1} - 2 y_n + y_{n+1}}{h^2} = -\lambda y_n \\
-\frac{2}{h^2}(y_0 - y_1) = -\lambda y_0 \\
y_N = y_{N-1}
\end{cases}
$$

## compile:
```bash
g++ main.cpp
```
## run:
```bash
./a.out 10
```