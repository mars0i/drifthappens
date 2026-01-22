(ns tips.fastmathMatrices
  (:require ;[clojure.math :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [fastmath.core :as fm]
            ;[fastmath.dev.codox :as codox]
            ;[fastmath.dev.clay :as utls]
            ;[fastmath.dev.plotly :as plt]
            ;[fastmath.dev.ggplot :as gg]
            ;[fastmath.stats :as stats]
           )
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))

;; Note that vector-only operations and matrix operations are
;; in different namespaces.

(def mat1 (mat/real-matrix [[1 2 3 4] [5 6 7 8]]))
(def mat2 (mat/mat         [[1 2 3 4] [5 6 7 8]]))

(def mat3 (mat/mat         [[1 2]
                            [3 4]]))

(def mat4 (mat/mat         [1 2 3 4]))

(def mat5 (mat/mat2x2 5))

(type mat1)
(type mat2)
(type mat3)
(type mat4)
(type mat5)

(mat/shape mat1)
(mat/shape mat2)
(mat/shape mat3)
(mat/shape mat4)
(mat/shape mat5)

(= mat1 mat2)

(Mat2x2. 1 2 3 4)
(seq (Mat2x2. 1 2 3 4))

(mat/entry mat1 1 1)
(mat/entry mat2 1 1)
(def mm1 (mat/mulmt mat1 mat2))
(mat/shape mm1)
(mat/entry mm1 1 1)

(clojure.pprint/pprint mm1)

(def mm2 (mat/mulm mat1 (mat/transpose mat2)))
(mat/shape mm2)
(mat/entry mm2 1 1)

(def mm3 (mat/mulm (mat/transpose mat1) mat2))
(mat/shape mm3)
(mat/entry mm3 1 1)
(mat/entry mm3 3 3)

(mat/mulm mat3 mat3)

(mat/mat [[1 2]])   ; 1x2 row matrix
(mat/mat [[1] [2]]) ; 2x1 column matrix
(mat/mat [1 2])     ; 2x1 column matric


(mat/mulm (mat/mat [[1 2 3] [4 5 6]]) (mat/mat [[10 100]]))

(mat/transpose (mat/mat [[1 2 3] [4 5 6]]))

;; ---

(def v2 (Vec2. -1 2.5))
;;(mat/shape v2)
(def v3 (Vec3. -1 2.5 17))
(def v6 (double-array (range 6)))

(v/mult v6 25) ; scalar mult of a vector
(v/emult v6 v6) ; element-wise, i.e. Hadamard mult


(mat/mulv (mat/mat [[1 2]
                    [3 4]])
          (double-array [10 100])) ; a double-array vector is treated as column

;; This fails btw:
;(mat/mulv (double-array [10 100])
;          (mat/mat [[1 2]
;                    [3 4]]))


(comment
  ;; Doesn't work.
  ;; How do you make a column vector? Is it a column matrix?
(mat/vtmul (double-array [10 100])
           (mat/mat [[1 2]
                     [3 4]]))

;; No
(mat/vtmul (mat/mat [[10 100]])
           (mat/mat [[1 2]
                     [3 4]]))

;; not this
(mat/vtmul (Vec2. 10 100)
           (mat/mat [[1 2]
                     [3 4]]))
)

;; OK so I don't know how to use vtmul, so the I can implement
;; mat/vector mult either with vectors on right, or by using 
;; matrices to represent vectors, in which case I can use either
;; direction.
