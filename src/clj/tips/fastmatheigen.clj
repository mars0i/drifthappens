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

(comment
  (fmat/shape small-tran-mat)
  (def smalldecomp (fmat/eigen-decomposition small-tran-mat))
  (def smalldecompc (fmat/eigen-decomposition small-tran-mat {:backend :colt}))
  (fmat/mat->array (fmat/decomposition-component smalldecomp :D))
  (fmat/mat->array (fmat/decomposition-component smalldecompc :D))
  (fmat/mat->array (fmat/decomposition-component smalldecomp :V))
  (fmat/mat->array (fmat/decomposition-component smalldecompc :V))
  (def evals(fmat/decomposition-component smalldecomp :real-eigenvalues))
  (def evalsc (fmat/decomposition-component smalldecompc :real-eigenvalues))
  (fmat/decomposition-component smalldecomp :imag-eigenvalues)
  (fmat/decomposition-component smalldecompc :imag-eigenvalues)

  (def evecs (fmat/decomposition-component smalldecomp :eigenvectors))
  (def evecsc (fmat/decomposition-component smalldecompc :eigenvectors))

  (fmat/mulv small-tran-mat (get evecsc 0))
  (fvec/mult (get evecsc 0) (get evalsc 0))

  (fmat/mulv small-tran-mat (get evecsc 1))
  (fvec/mult (get evecsc 1) (get evalsc 1))

  (fmat/mulv small-tran-mat (get evecsc 2))
  (fvec/mult (get evecsc 2) (get evalsc 2))


)

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

;pred-reprod-plots  ; display the plots


;; ---
; ### Plots from all three populations:

(comment
(def small-big-combo-plots 
  "States resulting from alternating plots from other sequences of plots."
  (interleave small-plots big-plots pred-reprod-plots))
)

;small-big-combo-plots   ; display the plots

