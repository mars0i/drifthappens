;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns scratch
  (:require ;[clojure.math :as m]
            ;[criterium.core :as crit]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]
            [fastmath.random :as frand]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.tableplot.v1.plotly :as plotly]
            [tablecloth.api :as tc]
            [drifthappens.wrightfisher :as wf]
            [utils.plotly :as uplot]
            [utils.fastmats :as mats]
            [utils.misc :as umisc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))

;; Example matrix from https://generateme.github.io/fastmath/clay/vector_matrix.html
(def M3x3 (Mat3x3. 1 2 3 -4 5 6 9 -8 7))

(def m3m (fmat/mulm M3x3 M3x3))

;; Generates a structure containing a map in which many components
;; are not delayed and not yet forced.  fmat/decomposition-component 
;; forces and replaces the :not-delivered value with the actual value.
(def m3m-edecomp (fmat/eigen-decomposition m3m))

(def v (fmat/decomposition-component m3m-edecomp :V))
(def d (fmat/decomposition-component m3m-edecomp :D))
(def vi (fmat/inverse v))

;; mulm takes only two matrix args
(reduce fmat/mulm [v d vi]) ; this product should be equal to m3m (up to float-slop)

v
;; A sequence of the columns of v:
(fmat/eigenvectors m3m)
(fmat/eigenvectors m3m true)
(map (fn [xs] (reduce + xs))
     (fmat/decomposition-component m3m-edecomp :eigenvectors))

d
(fmat/eigenvalues m3m) ; seq of pairs of [real, imag] components of eigenvals

;; Apparently these are not real and imaginary evals, but rather
;; real and imaginary components of eigenvalues:
(fmat/decomposition-component m3m-edecomp :real-eigenvalues) ; ? diag of d
(fmat/decomposition-component m3m-edecomp :imag-eigenvalues) ; ? sub/super pairs of d

;; Seems to be same as fmat/eigenvalues but as a matrix rather than
;; a clojure sequence.  Also seems to be same as :D from the
;; eigen-decomposition.
(fmat/eigenvalues-matrix m3m)

;; Seems to be same as :V in the eigen-decomposition:
(fmat/eigenvectors m3m)

;; --------------------------

;; From predreprod1.clj:

;; ---
;; ### Configuration and setup

(def fit-B 1.0)
(def big-fit-A 1.01)  ; large size has sel
(def small-fit-A 1.0) ; small size is pure drift
(def interval 16) ; interval between generations to display
(def half-interval (/ interval 2))

(def increments (iterate (partial + 16) 0))
(def generations increments)
(def half-generations (map (fn [n] (/ n 2)) increments))

(def num-gens 6)

;; These next two values should even numbers so that when divided,
;; we'll get integers:
(def small-N 2)
(def big-N 4)

(def half-small-N (/ small-N 2))
(def half-big-N (/ big-N 2))

;; ---
;; ### Small populations

(def small-pop-init 
  "A vector of initial probabilities of frequencies for a small population,
  typically with 1 as one element and zeros elswhere."
  (mats/mkvec (concat (repeat half-small-N 0.0) [1.0] (repeat half-small-N 0.0))))

(def small-tran-mat 
  "A single transition matrix for a small population."
  (wf/right-mult-tran-mat small-fit-A fit-B (dec (count small-pop-init))))

;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def small-tran-mats 
  "A sequence of M-to-M transition matrices, each of which is small-tran-mat raised to a power."
  (doall (mats/choose-mat-powers-separately small-tran-mat (take num-gens generations))))

(def small-prob-states 
  "States resulting from applying a product transition matrix to an initial state."
  (mats/make-prob-states small-tran-mats small-pop-init))

;; Plots made from the preceding sequence of states.
(def small-plots (mapv uplot/plot-both small-prob-states))

;small-plots  ; display the plots

;; ---
;; ### Large populations

(def big-pop-init 
  "A vector of initial probabilities of frequencies for a big population,
  typically with 1 as one element and zeros elswhere."
  (mats/mkvec (concat (repeat half-big-N 0.0) [1.0] (repeat half-big-N 0.0))))

(def big-tran-mat 
  "A single transition matrix for a big population."
  (wf/right-mult-tran-mat big-fit-A fit-B (dec (count big-pop-init)))) ; use fit-B for fit-A to make them equal

;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def big-tran-mats
  "A sequence of N-to-N transition matrices, each of which is big-tran-mat raised to a power."
  (doall (mats/choose-mat-powers-separately big-tran-mat (take num-gens generations))))

(def big-prob-states 
  "States resulting from applying a product transition matrix to an initial state."
  (mats/make-prob-states big-tran-mats big-pop-init))

;; Plots made from the preceding sequence of states.
(def big-plots (mapv uplot/plot-lines big-prob-states))

;big-plots  ; display the plots

;; ---
;; ### Alternating small/large populations

(def N (dec (count big-pop-init)))
(def M (dec (count small-pop-init)))

;; Shrinking and expanding matrices.
;;
;; Note M and N are swapped in the next two
(def predat-tran-mat
  "A transition matrix from a population of size N to a smaller population
  of size M resulting from predation."
  (wf/right-mult-tran-mat small-fit-A fit-B N M)) ; 

(def reprod-tran-mat
  "A transition matrix from a population of size M to a larger population
  of size N due to reproduction."
  (wf/right-mult-tran-mat big-fit-A fit-B M N)) ; use fit-B for fit-A to make them equal

;; Combine into one square matrix:
(def pred-reprod-mat 
  "A transition matrix from a populatiion of size N to a population of the
  same size. The matrix is the product of predat-tran-mat and reprod-drift
  mat.  It represents predation followed by reproductive growth."
  (fmat/mulm reprod-tran-mat predat-tran-mat))

;; The problem below with eigenvectors is this one:
;; https://github.com/generateme/fastmath/issues/42

(fmat/decomposition-component (fmat/eigen-decomposition m5x5) :ebal)

(comment
  ;; experiments with eigenvecs/vals
  pred-reprod-mat
  (type pred-reprod-mat)
  (type (fmat/array2d->RealMatrix pred-reprod-mat))
  (fmat/shape pred-reprod-mat)
  ;; I get "illegal state: convergence failed" errors, even for 100x100 matrix
  (fmat/eigenvalues pred-reprod-mat) ;; works with 5x5
  (fmat/eigenvectors pred-reprod-mat) ; fails.  Tried various conversions
  (def predreproddecomp (fmat/eigen-decomposition pred-reprod-mat))
  (fmat/decomposition-component predreproddecomp :D)
  (fmat/decomposition-component predreproddecomp :V)
  (type (fmat/decomposition-component predreproddecomp :V))
  (fmat/decomposition-component predreproddecomp :real-eigenvalues)
  (fmat/decomposition-component predreproddecomp :imag-eigenvalues)
  (fmat/decomposition-component predreproddecomp :eigenvectors) ; succeeds but empty seq

  ;; What happens if I convert to RealMat first?
  (def predreprodmat2d (fmat/array2d->RealMatrix pred-reprod-mat))
  (def ra2dcomp (fmat/eigen-decomposition predreprodmat2d))
  (fmat/eigenvectors predreprodmat2d) ; fails.  Tried various conversions
  (fmat/decomposition-component ra2dcomp :D)
  (fmat/decomposition-component ra2dcomp :V)
  (type (fmat/decomposition-component ra2dcomp :V))
  (fmat/decomposition-component ra2dcomp :real-eigenvalues)
  (fmat/decomposition-component ra2dcomp :imag-eigenvalues)
  (fmat/decomposition-component ra2dcomp :eigenvectors) ; succeeds but empty seq
)

(comment
;; We use half-generations here because each step involves two sampling processes.  
;; So each generation is analogous to two generations in the small and big models.
;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def pred-reprod-tran-mats
  "A sequence of N-to-N transition matrices, each of which is pred-reprod-mat raised to a power."
  (doall (mats/choose-mat-powers-separately pred-reprod-mat (take num-gens half-generations))))

(def pred-reprod-prob-states
  "States resulting from applying a product transition matrix to an initial state."
  (mats/make-prob-states pred-reprod-tran-mats big-pop-init))

;; Plots made from the preceding sequence of states.
(def pred-reprod-plots (mapv uplot/plot-both pred-reprod-prob-states))
)

(clojure.repl/dir fastmath.matrix)

;pred-reprod-plots  ; display the plots


;; ---
; ### Plots from all three populations:

(comment
(def small-big-combo-plots 
  "States resulting from alternating plots from other sequences of plots."
  (interleave small-plots big-plots pred-reprod-plots))
)

;small-big-combo-plots   ; display the plots


;; ----------------------------------------

;; double-double, toil and trouble
(def witch (mats/mkmat [[1 2 3 4]
                        [5 6 7 8]
                        [9 -1 -2 -3]
                        [-4 -5 -6 -7]]))

;; Note this takes four seqs, not a seq of seqs:
(def truck (fmat/rows->mat4x4 [1 2 3 4]
                              [5 6 7 8]
                              [9 -1 -2 -3]
                              [-4 -5 -6 -7]))

(def witchdecomp (fmat/eigen-decomposition witch))
(def truckdecomp (fmat/eigen-decomposition truck))

(def tv (fmat/decomposition-component witchdecomp :V))
(def wv (fmat/decomposition-component truckdecomp :V))
(def td (fmat/decomposition-component witchdecomp :D))
(def wd (fmat/decomposition-component truckdecomp :D))

(fmat/eigenvalues witch)
(fmat/eigenvalues truck)

(fmat/eigenvectors witch)
(fmat/eigenvectors truck)
(fmat/eigenvectors witch true)
(fmat/eigenvectors truck true)

;; ----------------------------------------


(def mdub (fm/seq->double-double-array [[1 2 3 4 5]
                                        [6 7 8 9 10]
                                        [-1 -2 -3 -4 -5]
                                        [-6 -7 -8 -9 -10]
                                        [1 2 3 4 5]]))
(type mdub)

(def m5x5 (fmat/rows->RealMatrix [[1 2 3 4 5]
                                   [6 7 8 9 10]
                                   [-1 -2 -3 -4 -5]
                                   [-6 -7 -8 -9 -10]
                                   [1 2 3 4 5]]))

(fmat/eigenvectors m5x5)
; (err) Wrong number of args (5) passed to: fastmath.matrix/rows->mat

(fmat/eigenvalues-matrix m5x5)
; (err) Wrong number of args (5) passed to: fastmath.matrix/rows->mat

(def m5x5-decomp (fmat/eigen-decomposition m5x5))

(fmat/decomposition-component m5x5-decomp :eigenvectors) ; => empty Clojure vector
(fmat/decomposition-component (fmat/eigen-decomposition m5x5) :eigenvectors) ; same thing

;; These succeed:
(fmat/eigenvalues m5x5)
(fmat/decomposition-component m5x5-decomp :real-eigenvalues)
(fmat/decomposition-component m5x5-decomp :imag-eigenvalues)

(fmat/decomposition-component m5x5-decomp :V)
; #object[org.apache.commons.math3.linear.Array2DRowRealMatrix 0xc83fff9 "Array2DRowRealMatrix
; {{0.5379033159,-0.1179138184,1.4501010624,-1.0040216581,3.293351034},
;  {0.4882862549,-0.9871935558,-0.8123428453,1.7903179275,-7.7549758304},
;  {-0.5379033159,0.1179138184,-1.3828351958,-0.2704195191,3.4573201361},
;  {-0.4882862549,0.9871935558,-0.5977053223,-0.8140281116,3.1768830828},
;  {0.5379033159,-0.1179138184,1.3427823009,0.2981513614,-2.1725784226}}"]

(fmat/decomposition-component m5x5-decomp :D)
;#object[org.apache.commons.math3.linear.Array2DRowRealMatrix 0x276a9d93 "Array2DRowRealMatrix
;{{0.5,3.1224989992,0.0,0.0,0.0},
; {-3.1224989992,0.5,0.0,0.0,0.0},
; {0.0,0.0,0.0,0.0,0.0},
; {0.0,0.0,0.0,0.0,0.0},
; {0.0,0.0,0.0,0.0,0.0}}"]

