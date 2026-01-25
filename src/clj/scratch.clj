;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns scratch
  (:require ;[clojure.math :as m]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]
            ;[fastmath.dev.codox :as codox]
            ;[fastmath.dev.clay :as utls]
            ;[fastmath.dev.plotly :as plt]
            ;[fastmath.dev.ggplot :as gg]
            ;[fastmath.stats :as stats]
            [scicloj.kindly.v4.kind :as kind]
            [utils.misc :as misc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))


(def v2 (Vec2. -1 2.5))
(type v2)
(def v3 (Vec3. -1 2.5 -3.25))
(type v3)
(def v3' (fvec/vec3 -1 2.5 -3.25)) ; same thing as preceding
(type v3')
(def v4 (Vec4. -1 0.0 -3.25 4))
(type v4)
(def av5 (fvec/array-vec [-1 0.0 -3.25 4 -23]))
(type av5)
(def av4 (fvec/make-vector 4 [-1 0.0 -3.25 4 -23])) ; ignores extra elements, uses 0.0 for unprovided elements
(type av4)
(def av4 (fvec/make-vector 7 [-1 0.0 -3.25 4 -23])) ; ignores extra elements, uses 0.0 for unprovided elements
fvec/generate-vec3
fvec/array->vec3
fvec/seq->vec3
make-array
(def da10 (double-array (range 10))) ; makes a Java array of doubles (?)
(type da10)


(def mat1 (fmat/real-matrix [[1 2 3 4]
                             [5 6 7 8]]))

(def mat2 (fmat/mat         [[1 2 3]
                             [4 5 6]]))

(def mat3 (fmat/mat         [[1 2]
                             [5 6]]))
(type mat1)
(type mat2)
(type mat3)

(def colmat3 (fmat/mat [[10]
                       [100]
                       [1000]]))

(def rowmat3 (fmat/mat [[10 100 1000]]))

(def rowmat2 (fmat/mat [[10 100]]))

(def da3 (double-array [10 100 1000]))
(type da3)
(def la3 (long-array [10 100 1000]))
(type la3)

(def v3 (fvec/vec3 10 100 1000))

(fmat/mulm mat2 colmat3)
(fmat/mulm rowmat2 mat2)

;(fmat/mulv mat2 colmat3) ; fails
;(fmat/mulv mat2 rowmat3) ; fails
(fmat/mulv mat2 da3) ; succeeds
;(fmat/mulv mat2 la3) ; fails
;(fmat/mulv mat2 v3) ; fails
;; So it seems that mulv only works when the vector is represented by a Java double array.


;; Definitions from https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(def M2x2 (Mat2x2. 1 2 3 4))
(type M2x2)
(def M3x3 (Mat3x3. 1 2 3 -4 5 6 9 -8 7))
(type M3x3)
(def M4x4 (Mat4x4. 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16))
(type M4x4)
(def RealMat (fmat/real-matrix [[1 2 3] [4 5 6]]))
(type RealMat)

;; Note that row representations are treated as column vectors here:
(fmat/mulv M2x2 (fvec/vec2 10 20))
(type (fmat/mulv M2x2 (fvec/vec2 10 20)))
(fmat/vtmul M2x2 (fvec/vec2 10 20))
;; NOTE PRECEDING: It reps mult on the left by the vec, on right by mat, 
;; but the vec arg is on the right, and is treated as a col vec that will
;; be transposed into a row vec.

;(fmat/mulv RealMat (fvec/vec3 1 2 3)) ; fails
(fmat/mulv RealMat (fvec/vec->RealVector (fvec/vec3 1 2 3)))
(fmat/mulv RealMat (fvec/vec->RealVector (range 1 4)))
(type (fmat/mulv RealMat (fvec/vec->RealVector (fvec/vec3 1 2 3))))
(type (fvec/vec->RealVector (fvec/vec3 1 2 3)))
(type (fmat/mulv RealMat (fvec/vec->RealVector (range 1 4))))

(fmat/vtmul RealMat (fvec/vec->RealVector [1 2])) 

(def rv2 (fvec/vec->RealVector [10 100])) ; not seq-able
(def rv3 (fvec/vec->RealVector [10 100 1000])) ; not seq-able
(def da2 (double-array [10 100])) ; seq-able
(def da3 (double-array [10 100 1000])) ; seq-able

(def av2 (fvec/array-vec [10 100])) ; seq-able
(def av3 (fvec/array-vec [10 100 1000])) ; seq-able
(type av3)

(fmat/mulv mat2 rv3)   ; works correctly
(fmat/mulv mat2 da3)   ; works correctly
;(fmat/mulv mat2 av3)   ; fails
(fmat/vtmul mat2 rv2)  ; works correctly. See note above about semantics of vtmul.
(fmat/vtmul mat2 da2)  ; works correctly.
;(fmat/vtmul mat2 av2)  ; fails
;; SUMMARY:
;; Fastmath vector/matrix mult works only with apache RealVector and Java double arrays.
;; It doesn't work with other arbitrary-length fastmath vector
;; representations such as fastmath.vector.ArrayVec (in fastmath 3.0.0-alpha, at least)

;; btw:
;; fvec/seq->RealVector ;; doesn't exist
(def rv3' (fvec/vec->RealVector (list 10 100 1000))) ; succeeds but I think the result is wrong
(def rv3'' (fvec/vec->RealVector (map identity (list 10 100 1000)))) ; succeeds but I think the result is wrong
