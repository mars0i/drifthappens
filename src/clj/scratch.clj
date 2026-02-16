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
(fmat/decomposition-component m3m-edecomp :eigenvectors)

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

