# Class 6: R Functions
Katherine Quach (A18541014)

- [Background](#background)
- [Our First Function](#our-first-function)
- [A Second Function](#a-second-function)
- [A Protein generating function](#a-protein-generating-function)

## Background

All functions in R have at least 3 things:

- A **name** that we use to call the function
- One or more input **arguments**
- The **body** the lines of R code that do the work

## Our First Function

Let’s write a silly wee function called `add()` to add some numbers (the
input arguments)

``` r
add <-  function(x, y){
  x + y
}
```

Now we can use this function

``` r
add(x = 100, y = 1)
```

    [1] 101

``` r
add(x = 10, y = 10)
```

    [1] 20

``` r
add(x = c(100, 1, 100), y = 1)
```

    [1] 101   2 101

> Q. What if I give a multiple element vector to `x` and `y`

``` r
add(x = c(100, 1), y = c(100, 1))
```

    [1] 200   2

> Q. What if I give three inputs to the function

``` r
#add(x = c(100, 1), y = 1, z = 1)
```

> Q. What if I give only one input to the add function

``` r
addnew <-  function(x, y = 1){
  x + y
}
```

``` r
addnew(x = 100)
```

    [1] 101

``` r
addnew(c(100, 1), 100)
```

    [1] 200 101

If we write our function with input arguments having no default value
then the user will be required to set them when they use the function.
We can give our input arguments “default” values by setting them equal
to some sensible value - e.g. y = 1 in the `addnew()` function

## A Second Function

Let’s try something more interesting: Make a sequence generating tool…

The `sample()` function can be a useful starting point here:

``` r
sample(1:10, size = 4)
```

    [1]  3 10  9  8

> Q. Generate 9 random numbers taken from the input vector x = 1:10?

``` r
sample(1:10, size = 9)
```

    [1]  8  7  1 10  9  2  6  3  5

> Q. Generate 12 random numbers taken from the input vector x = 1:10?

``` r
sample(1:10, size = 12, replace = TRUE)
```

     [1]  3 10 10  5  8  2  4  7  6  3 10  5

> Q. Write code for the `sample()` function that generates nucleotide
> sequences of length 6

``` r
sample(c("A", "G", "T", "C"), size = 6, replace = TRUE)
```

    [1] "A" "G" "A" "C" "C" "A"

> Q. Write a first function `generate_dna()` that returns a *user
> specified length* DNA sequence:

``` r
generate_dna <- function(len = 6){
  sample(c("A", "G", "T", "C"), size = len, replace = TRUE)
}
```

``` r
generate_dna(len = 100)
```

      [1] "C" "G" "T" "C" "A" "T" "A" "C" "T" "G" "T" "A" "A" "G" "T" "A" "T" "C"
     [19] "G" "G" "T" "A" "C" "A" "A" "C" "C" "A" "A" "T" "G" "G" "G" "C" "A" "G"
     [37] "C" "T" "C" "G" "A" "G" "C" "A" "G" "C" "C" "C" "G" "T" "A" "G" "A" "C"
     [55] "C" "C" "C" "T" "T" "A" "T" "C" "C" "G" "T" "T" "T" "A" "C" "A" "C" "T"
     [73] "C" "G" "T" "G" "C" "C" "C" "T" "T" "T" "A" "G" "C" "C" "A" "A" "G" "G"
     [91] "T" "A" "A" "C" "C" "A" "A" "T" "A" "A"

> **Key Points** Every function in R looks fundamentally the same in
> terms of its structure. Basically 3 things: name, input, and body

    name <- function(input) {
    body
    }

> Functions can have multiple inputs. These can be **required**
> arguments or **optional** arguments. With optional arguments having a
> set of default value

> Q. Modify and improve our `generate_dna()` function to return its
> generated sequence in a more standard format like “AGTAGTA” rather
> than the vector “A”, “C”, “G”, “A”

``` r
generate_dna <- function(len = 6, fasta = TRUE){
  ans <- sample(c("A", "G", "T", "C"), 
         size = len, replace = TRUE)
  if(fasta) {
    cat("Single-element vector output")
    ans <- paste(ans, collapse = "")
  } else {
    cat("Multi-element vector output")
  }
return(ans)
  }
generate_dna(fasta = FALSE)
```

    Multi-element vector output

    [1] "G" "T" "T" "T" "G" "C"

``` r
generate_dna(fasta = TRUE)
```

    Single-element vector output

    [1] "TTTAGA"

The `paste()` function - its job is to join up or stick together (a.k.a.
paste) input strings together

``` r
paste("alice", "loves R", sep = " ")
```

    [1] "alice loves R"

Flow Control means where the R brain goes in your code

``` r
good_mood <- FALSE

if(good_mood) {
  cat("Great")
} else {
  cat("Bummer!")
}
```

    Bummer!

## A Protein generating function

> Q. Write a function, called `generate_protein()` that generates a user
> specified length protein sequence

There are 20 natural amino-acids

``` r
aa <- c("A", "R", "N", "D", "C", "Q", "E",
        "G", "H", "I", "L", "K", "M", "F",
        "P", "S", "T", "W", "Y", "V")
```

``` r
generate_protein <- function(len) {
  # The amino-acids to sample from
   aa <- c("A", "R", "N", "D", "C", "Q", "E",
        "G", "H", "I", "L", "K", "M", "F",
        "P", "S", "T", "W", "Y", "V")
   # Draw n = len amino-acids to make our sequence
  ans <- sample(aa, size = len, replace = T)
  ans <- paste(ans, collapse = "")
  return(ans)
}
```

``` r
myseq <- generate_protein(42)
myseq
```

    [1] "VDMEYHNYVCIVTKELHGISRDPAAQYQNHPSHDTDEEKGRC"

> Q. Use that function to generate random protein sequences between
> length 6 and 12

``` r
generate_protein(6)
```

    [1] "KGGGTH"

``` r
generate_protein(7)
```

    [1] "MNDWPWH"

``` r
generate_protein(8)
```

    [1] "CPDKTIPH"

``` r
generate_protein(9)
```

    [1] "PGLKAKLFF"

``` r
generate_protein(10)
```

    [1] "GEVKEQPYGI"

``` r
generate_protein(11)
```

    [1] "SDVSKVCWKNS"

``` r
generate_protein(12)
```

    [1] "YACYHEKPLSME"

``` r
for(i in 6:12) {
  # FASTA ID line ">id"
  cat(">", i, "\n")
  # Protein sequence line
  cat(generate_protein(i), "\n")
}
```

    > 6 
    FWNPRK 
    > 7 
    QQDNDGM 
    > 8 
    RTMNSFEY 
    > 9 
    PGPMWCDCT 
    > 10 
    MTNMTIRNGY 
    > 11 
    KCGRSTFLNYP 
    > 12 
    IKDWFNAGLQCP 

> Q. Are any of your sequences unique i.e. not found anywhere in nature?

These are all of the unique sequences that are not found anywhere in
nature \> 8 WFQVVCQA \> 9 RWMACNRKH \> 10 PNFQKMMKCK \> 11 DHPGKTFTIFA
\> 12 GHMQEYEIGLAR
