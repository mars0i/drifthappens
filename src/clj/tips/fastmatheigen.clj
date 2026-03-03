;; Illustrations of extraction of eigenvalues and eigenvectors
;; using fastmath on my drift matrices.
(ns fastmatheigen
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

;; --------------------------

;; From predreprod1.clj:

;; ---
;; ### Configuration and setup

(def fit-B 1.0)
(def big-fit-A 1.0)  ; large size has sel
(def small-fit-A 1.0) ; small size is pure drift
(def interval 16) ; interval between generations to display
(def half-interval (/ interval 2))

;(def increments (iterate (partial + 16) 0))
(def increments (iterate (partial + 1) 0))
(def generations increments)
(def half-generations (map (fn [n] (/ n 2)) increments))

(def num-gens 100)

;; These next two values should even numbers so that when divided,
;; we'll get integers:
(def small-N 10)
(def big-N 100)

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
(def small-plots (mapv uplot/plot-dots-lines small-prob-states))

;small-plots  ; display the plots

;; ---
;; ### Large populations

;(def big-pop-init 
;  "A vector of initial probabilities of frequencies for a big population,
;  typically with 1 as one element and zeros elswhere."
;  (mats/mkvec (concat (repeat half-big-N 0.0) [1.0] (repeat half-big-N 0.0))))

;; alt init cond doesn't seem to work
;(def big-pop-init (mats/mkvec (assoc (vec (repeat (inc big-N) 0)) 5 1)))

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


;; We use half-generations here because each step involves two sampling processes.  
;; so each generation is analogous to two generations in the small and big models.
;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def pred-reprod-tran-mats
  "A sequence of N-to-N transition matrices, each of which is pred-reprod-mat raised to a power."
  (doall (mats/choose-mat-powers-separately pred-reprod-mat (take num-gens half-generations))))

(def pred-reprod-prob-states
  "States resulting from applying a product transition matrix to an initial state."
  (mats/make-prob-states pred-reprod-tran-mats big-pop-init))

;; Plots made from the preceding sequence of states.
(def pred-reprod-plots (mapv uplot/plot-dots-lines pred-reprod-prob-states))

;pred-reprod-plots  ; display the plots


;; ---
; ### Plots from all three populations:

(def small-big-combo-plots 
  "States resulting from alternating plots from other sequences of plots."
  (interleave small-plots big-plots pred-reprod-plots))

small-big-combo-plots   ; display the plots

;; -------------------------------
;; ### Experiments

(comment
  (fmat/shape small-tran-mat)

  ;; Apache decomposition:
  (def smalldecomp (fmat/eigen-decomposition small-tran-mat))
  (fmat/mat->array (fmat/decomposition-component smalldecomp :D))
  (fmat/mat->array (fmat/decomposition-component smalldecomp :V))
  (def evals (fmat/decomposition-component smalldecomp :real-eigenvalues))
  (fmat/decomposition-component smalldecomp :imag-eigenvalues)
  (def evecs (fmat/decomposition-component smalldecomp :eigenvectors))
  (map fvec/vec->seq evecs)

  ;; A single multiplication of the initial state by the matrix:
  (fmat/mulv small-tran-mat small-pop-init)
  ;; Now can I do the same thing with the deomposition?
  ;; First I need to rep the vec (small-pop-init) in terms of 
  ;; the eigenvectors.


  ;; Colt decomposition:
  (def smalldecompc (fmat/eigen-decomposition small-tran-mat {:backend :colt}))
  (fmat/mat->array (fmat/decomposition-component smalldecompc :D))
  (fmat/mat->array (fmat/decomposition-component smalldecompc :V))

  (fmat/decomposition-component smalldecompc :imag-eigenvalues)
  (def evalsc (fmat/decomposition-component smalldecompc :real-eigenvalues))
  ;; Ewens p. 22 says that the leading non-unit eigenvalue is
  ;; $1-1/(2N)$ for a diploid model  This should be $1-1/N$ for
  ;; haploid, or $(N-1)/N$.
  (double (- 1 (/ small-N)))
  ;; So for $N=10$, this is 0.9.  With $N=10$, the first two
  ;; eigenvalues from both Apache and Colt above are approx 1,
  ;; and the third eigenvalue is indeed approx 0.9.

  ;; The first eigenvector is (1,0,...,0).  The second consists of
  ;; negative numbers, with -0.04 in the first place, very small negative
  ;; numbers in most of the others, and -1 in the last (mod float slop).
  (def evecsc (fmat/decomposition-component smalldecompc :eigenvectors))
  ;; The first eigenvector is (1,0,...,0).  The second consists of
  ;; negative numbers, with -0.49 in the first place, very small negative
  ;; numbers in most of the others, and -1 in the last (mod float slop).
)

(comment
  (fmat/shape big-tran-mat)

  ;; Apache decomposition doesn't work for this matrix:
  (def bigdecomp (fmat/eigen-decomposition big-tran-mat))

  ;; Colt decomposition:
  (def bigdecompc (fmat/eigen-decomposition big-tran-mat {:backend :colt}))
  (fmat/mat->array (fmat/decomposition-component bigdecompc :D))
  (fmat/mat->array (fmat/decomposition-component bigdecompc :V))

  ;; Note all of the imaginary components are zero, but the nonzero ones 
  ;; are extremely small, so I'll ignore them:
  (fmat/decomposition-component bigdecompc :imag-eigenvalues)
  (def evalsc (fmat/decomposition-component bigdecompc :real-eigenvalues))

  ;; Ewens p. 22 says that the leading non-unit eigenvalue is
  ;; $1-1/(2N)$ for a diploid model  This should be $1-1/N$ for
  ;; haploid, or $(N-1)/N$.
  ;; So for $N=100$, this is 0.99.  With $N=100$, the first eigenvalue
  ;; is exactly 1.0.  The second is close: 0.9999999999999993. 
  ;; The third eigenvalue is 0.9899999999999993, which does round to 0.99.

  (def evecsc (fmat/decomposition-component bigdecompc :eigenvectors))
  (evecsc 0) ; 1 followed by 100 zeros.
  (evecsc 1) ; last elem approx 1, first is 6.785, middles are very small.
  (evecsc 2) ; first and last elem are approx 4.97; others are near -0.1, or -0.08 at either end
)

(comment
  ;; experiments with eigenvecs/vals
  pred-reprod-mat
  (type pred-reprod-mat)
  (fmat/shape pred-reprod-mat)

  ;; By default these functions use Apache's math3.linear.SchurTransformer/transform
  ;; This fails for some of my matrices.  Generateme added other back ends.  
  ;; Adding {:backend :colt} to `eigen-decomposition` works, and Generateme
  ;; says they'll add this to `eigenvalues` and `eigenvectors`.
  ;(fmat/eigenvalues pred-reprod-mat) ;; works with 5x5
  ;(fmat/eigenvectors pred-reprod-mat) ; fails.  Tried various conversions
  ;(def predreproddecomp (fmat/eigen-decomposition pred-reprod-mat)) ; 

  (def predreproddecomp (fmat/eigen-decomposition pred-reprod-mat {:backend :colt}))

  (fmat/mat->array 
    (fmat/decomposition-component predreproddecomp :D)
  )
  (fmat/mat->array 
    (fmat/decomposition-component predreproddecomp :V)
  )
  (type (fmat/decomposition-component predreproddecomp :V))
  (fmat/decomposition-component predreproddecomp :real-eigenvalues)
  (fmat/decomposition-component predreproddecomp :imag-eigenvalues)
  (def evals (fmat/decomposition-component predreproddecomp :eigenvectors))
  (map type evals)
  (map fvec/vec->seq evals)
  (count evals)
  (evals 0)
  (get (evals 100) 100)
  (map #(get % 100) evals)

  ;; ---

  ;; What happens if I convert to RealMat first?
  (def pred-reprod-2d (fmat/array2d->RealMatrix pred-reprod-mat))
  (type pred-reprod-2d)
  (fmat/shape pred-reprod-2d)
  (fmat/mat->array2d pred-reprod-2d)

  (def ra2dcomp (fmat/eigen-decomposition pred-reprod-2d))
  (def ra2dcomp (fmat/eigen-decomposition pred-reprod-2d {:backend :colt}))
  (def otherdecomp (fmat/eigen-decomposition (fmat/array2d->RealMatrix big-tran-mat) {:backend :colt})) ;; ALTERNATIVE TEST
  (def evals (fmat/decomposition-component ra2dcomp :real-eigenvalues))
  (def ievals (fmat/decomposition-component ra2dcomp :imag-eigenvalues))
  (def evecs (fmat/decomposition-component ra2dcomp :eigenvectors))
  (fmat/eigenvectors pred-reprod-2d) ; fails.  Tried various conversions
  (fmat/eigenvalues pred-reprod-2d) ; fails.  Tried various conversions
  (fmat/decomposition-component ra2dcomp :D)
  (fmat/decomposition-component ra2dcomp :V)
  (type (fmat/decomposition-component ra2dcomp :V))
  (type evals) ; evals is Java vector not a Clojure vector
  (map fvec/vec->seq evecs)

  ;; test the eigenvector/value relationship for each pair:
  ;; Only works with real eigenvalues
  (map #(fvec/delta-eq (fmat/mulv pred-reprod-2d %1) ; matrix-vector mult
                       (fvec/mult %1 %2)) ; vector-scalar mult
       evecs
       (vec evals))
)

