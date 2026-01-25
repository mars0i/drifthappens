(ns mwe
  (:require [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat])
  (:import [fastmath.vector Vec2 Vec3]
           [fastmath.matrix Mat2x2 Mat3x3]))


;; fmat/real-matrix does the same thing:
(def m (fmat/mat [[1 2 3]
                  [4 5 6]]))

;; All of these four vector reps are seq-able, except rv3. 
;; So double-array is the only rep that is both seq-able and works with
;; matrix/vector multiplication.
(def da3 (double-array [1000 100 10])) ; double/1
(def rv3 (fvec/vec->RealVector [1000 100 10])) ; org.apache.commons.math3.linear.ArrayRealVector
(def v3 (fvec/vec3 1000 100 10)) ; fastmath.vector.Vec3
(def av3 (fvec/array-vec [1000 100 10])) ; fastmath.vector.ArrayVec


(comment
  (fmat/mulv m da3) ; works
  (fmat/mulv m rv3) ; works

  (fmat/mulv m v3)  ; err: class fastmath.vector.Vec3 cannot be cast to class [D 
  ;; But this fastmath.vector.Vec3 can be multiplied by a fastmath.matrix.Mat3x3

  (fmat/mulv m av3) ; err: class fastmath.vector.ArrayVec cannot be cast to class [D

  ;; The error messages also say things like:
  ;; (fastmath.vector.Vec3 is in unnamed module of loader 
  ;; clojure.lang.DynamicClassLoader @53a707a2; [D is in
  ;;  module java.base of loader 'bootstrap')
)

(def da2 (double-array [100 10]))
(def rv2 (fvec/vec->RealVector [100 10]))
(def v2 (fvec/vec2 100 10))
(def av2 (fvec/array-vec [100 10]))

(def m2x2 (fmat/mat2x2 1 2
                       3 4))

(comment
  (fmat/vtmul m da2) ; works correctly.
  (fmat/vtmul m rv2) ; works correctly.
  (fmat/vtmul m v2)  ; fails
  (fmat/vtmul m av2) ; fails

  (fmat/vtmul m2x2 v2)  ; works correctly.
  (fmat/vtmul m2x2 av2) ; fails
)

(comment
  (type da2)
  (type rv2)
  (type v2)
  (type av2)

  (type da3)
  (type rv3)
  (type v3)
  (type av3)
)

(def m3x3 (fmat/mat3x3 1 2 3
                       4 5 6
                       7 8 9))

;; SUMMARY:
;; Fastmath vector/matrix mult works only with apache RealVector and Java double arrays.
;; It doesn't work with other arbitrary-length fastmath vector
;; representations such as fastmath.vector.ArrayVec (in fastmath 3.0.0-alpha, at least)
