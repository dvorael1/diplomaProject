library('FunChisq')


# Not run: 
# In all examples, x is the row variable and y is the column
#    variable of a table.

# Example 1. Simulating a noisy function where y=f(x),
#            x may or may not be g(y)

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="functional",
                        noise=0.2, n.tables = 1,
                        row.marginal = c(0.3,0.2,0.3,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 1. Functional pattern")
plot_table(tbls$sample.list[[1]], main="Ex 1. Sampled pattern (noise free)")
plot_table(tbls$noise.list[[1]], main="Ex 1. Sampled pattern with 0.2 noise")
plot.new()

# Example 2. Simulating a noisy functional pattern where
#            y=f(x), x may or may not be g(y)

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="functional",
                        noise=0.5, n.tables = 1,
                        row.marginal = c(0.3,0.2,0.3,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 2. Functioal pattern", col="seagreen2")
plot_table(tbls$sample.list[[1]], main="Ex 2. Sampled pattern (noise free)", col="seagreen2")
plot_table(tbls$noise.list[[1]], main="Ex 2. Sampled pattern with 0.5 noise", col="seagreen2")
plot.new()


# Example 3. Simulating a noise free many.to.one function where
#            y=f(x), x!=f(y).

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="many.to.one",
                        noise=0.2, n.tables = 1,
                        row.marginal = c(0.4,0.3,0.1,0.2))
par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 3. Many-to-one pattern", col="limegreen")
plot_table(tbls$sample.list[[1]], main="Ex 3. Sampled pattern (noise free)", col="limegreen")
plot_table(tbls$noise.list[[1]], main="Ex 3. Sampled pattern with 0.2 noise", col="limegreen")
plot.new()

# Example 4. Simulating noise-free discontinuous
#   pattern where y=f(x), x may or may not be g(y)

tbls <- simulate_tables(n=100, nrow=4, ncol=5,
                        type="discontinuous", noise=0.2,
                        n.tables = 1, row.marginal = c(0.2,0.4,0.2,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 4. Discontinuous pattern", col="springgreen3")
plot_table(tbls$sample.list[[1]], main="Ex 4. Sampled pattern (noise free)", col="springgreen3")
plot_table(tbls$noise.list[[1]], main="Ex 4. Sampled pattern with 0.2 noise", col="springgreen3")
plot.new()


# Example 5. Simulating noise-free dependent.non.functional
#            pattern where y!=f(x) and x and y are statistically
#            dependent.

tbls <- simulate_tables(n=100, nrow=4, ncol=5,
                        type="dependent.non.functional", noise=0.3,
                        n.tables = 1, row.marginal = c(0.2,0.4,0.2,0.2))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 5. Dependent.non.functional pattern",
           col="sienna2", highlight="none")
plot_table(tbls$sample.list[[1]], main="Ex 5. Sampled pattern (noise free)",
           col="sienna2", highlight="none")
plot_table(tbls$noise.list[[1]], main="Ex 5. Sampled pattern with 0.3 noise",
           col="sienna2", highlight="none")
plot.new()

# Example 6. Simulating a pattern where x and y are
#            statistically independent.

tbls <- simulate_tables(n=100, nrow=4, ncol=5, type="independent",
                        noise=0.3, n.tables = 1,
                        row.marginal = c(0.4,0.3,0.1,0.2),
                        col.marginal = c(0.1,0.2,0.4,0.2,0.1))

par(mfrow=c(2,2))
plot_table(tbls$pattern.list[[1]], main="Ex 6. Independent pattern",
           col="cornflowerblue", highlight="none")
plot_table(tbls$sample.list[[1]], main="Ex 6. Sampled pattern (noise free)",
           col="cornflowerblue", highlight="none")
plot_table(tbls$noise.list[[1]], main="Ex 6. Sampled pattern with 0.3 noise",
           col="cornflowerblue", highlight="none")
plot.new()


## Not run: 
# Example 1. Asymptotic functional chi-square test
x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
fun.chisq.test(x) # strong functional dependency
fun.chisq.test(t(x)) # weak functional dependency

# Example 2. Normalized functional chi-square test
x <- matrix(c(8,0,8,0,8,0,2,0,2), 3)
fun.chisq.test(x, method="nfchisq") # strong functional dependency
fun.chisq.test(t(x), method="nfchisq") # weak functional dependency

# Example 3. Exact functional chi-square test
x <- matrix(c(4,0,4,0,4,0,1,0,1), 3)
fun.chisq.test(x, method="exact") # strong functional dependency
fun.chisq.test(t(x), method="exact") # weak functional dependency

# Example 4. Exact functional chi-square test on a real data set
#            (Shen et al., 2002)
# x is a contingency table with row variable for p53 mutation and
#   column variable for CIMP
x <- matrix(c(12,26,18,0,8,12), nrow=2, ncol=3, byrow=TRUE)

# Test the functional dependency: p53 mutation -> CIMP
fun.chisq.test(x, method="exact")

# Test the functional dependency CIMP -> p53 mutation
fun.chisq.test(t(x), method="exact")

# Example 5. Asymptotic functional chi-square test with simulated distribution
x <- matrix(c(20,0,20,0,20,0,5,0,5), 3)
fun.chisq.test(x, method="simulate.p.value")
fun.chisq.test(x, method="simulate.p.value", simulate.n = 1000)

## End(Not run)

## End(Not run)
