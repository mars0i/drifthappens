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


(def A (mats/mkmat [[0 5 -10 ][ 0 22 16 ][ 0 -9 -2]]))
(fmat/mulv A (mats/mkvec [-5 -4 3]))
(fmat/mulv A (mats/mkvec [1 0 0]))

(def fow (fmat/eigen-decomposition A))
(:real-eigenvalues fow)

(def v (fmat/decomposition-component fow :V))
(def d (fmat/decomposition-component fow :D))
(def di (fmat/inverse d)) ; noninvertible ; nonindvertible
;(fmat/mulm v d (fmat/inverse d))

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

;; A sequence of the columns of v:
(fmat/decomposition-component m3m-edecomp :eigenvectors)
v

(fmat/decomposition-component m3m-edecomp :real-eigenvalues) ; ? diag of d
(fmat/decomposition-component m3m-edecomp :imag-eigenvalues) ; ? sub/super pairs of d
d

(fmat/decomposition-component m3m-edecomp :complex?)

;; see also https://generateme.github.io/fastmath/clay/vector_matrix.html#eigenvalues-and-singular-values
