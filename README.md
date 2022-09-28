# eigenvalue-problem

## problem:

$$
\Huge{}
\begin{cases}
\frac{y_{n+1} - 2 \, y_n + y_{n+1}}{h^2} = -\lambda y_n \\
-2 \, \frac{y_0 - y_1}{h^2} = -\lambda y_0 \\
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