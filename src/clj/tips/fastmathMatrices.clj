(ns tips.fastmathMatrices
  (:require ;[clojure.math :as m]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]
            ;[fastmath.dev.codox :as codox]
            ;[fastmath.dev.clay :as utls]
            ;[fastmath.dev.plotly :as plt]
            ;[fastmath.dev.ggplot :as gg]
            ;[fastmath.stats :as stats]
           )
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))



;; At https://clojurians.zulipchat.com/#narrow/channel/151924-data-science/topic/Matrix-vector.20multiplication.20in.20fastmath/near/569973078
;; @generateme recommends using `fastmath.core/seq->double-double-array`
;; to create matrices, and `fastmath.core/sec->double-array` to create
;; (column) vectors. 
;; i.e. Use `seq->double-double-array` to create (Java) matrices rather
;; than fastmath.matrix/mat, for example, and use `seq->double-array` or
;; Clojure's `double-array` to create vectors rather than other methods.
;; This will allow fastmath.matrix/multv, vtmult, and multm to work, and
;; double arrays are seq-able.
;; See also: https://generateme.github.io/fastmath/clay/core.html#array-and-sequence-conversion

;; ---------------------------------------

(def da1 (double-array [1 2 3]))
(def da2 (fm/seq->double-array [1 2 3]))

(def da1copy (double-array da1))
(def da1same (fm/seq->double-array da1))

(= da1 da2)
(= da1 da1copy)
(= da1 da1same)
(= (seq da1) (seq da2))

(def dam (fm/seq->double-double-array [[1 2 3]
                                       [4 5 6]
                                       [7 8 9]
                                       [10 11 12]]))
(type dam)

(def mm (fmat/mat [[1 2 3]
                   [4 5 6]
                   [7 8 9]
                   [10 11 12]]))
(type mm)

;; ---------------------------------------

(def m23 (fm/seq->double-double-array [[1 2 3]
                                       [4 5 6]]))
(def v2 (fm/seq->double-array [100 10]))
(def v3 (double-array [1000 100 10]))

;; ---------------------------------------

;; Note that vector-only operations and matrix operations are
;; in different namespaces.

(def mat1 (fmat/real-matrix [[1 2 3 4] [5 6 7 8]]))
(def mat2 (fmat/mat         [[1 2 3 4] [5 6 7 8]]))

(def mat3 (fmat/mat         [[1 2]
                            [3 4]]))

(def mat4 (fmat/mat         [1 2 3 4]))

(def mat5 (fmat/mat2x2 5))

(type mat1)
(type mat2)
(type mat3)
(type mat4)
(type mat5)

(fmat/shape mat1)
(fmat/shape mat2)
(fmat/shape mat3)
(fmat/shape mat4)
(fmat/shape mat5)

(= mat1 mat2)

(Mat2x2. 1 2 3 4)
(seq (Mat2x2. 1 2 3 4))

(fmat/entry mat1 1 1)
(fmat/entry mat2 1 1)
;(def mm1 (fmat/mulmt mat1 mat2))
;(fmat/shape mm1)
;(fmat/entry mm1 1 1)

;(clojure.pprint/pprint mm1)

;(def mm2 (fmat/mulm mat1 (fmat/transpose mat2)))
;(fmat/shape mm2)
;(fmat/entry mm2 1 1)

;(def mm3 (fmat/mulm (fmat/transpose mat1) mat2))
;(fmat/shape mm3)
;(fmat/entry mm3 1 1)
;(fmat/entry mm3 3 3)

;(fmat/mulm mat3 mat3)

(fmat/mat [[1 2]])   ; 1x2 row matrix
(fmat/mat [[1] [2]]) ; 2x1 column matrix
(fmat/mat [1 2])     ; 2x1 column matric


(fmat/transpose (fmat/mat [[1 2 3] [4 5 6]]))

;; ---

(def v2 (Vec2. -1 2.5))
;;(fmat/shape v2)
(def v3 (Vec3. -1 2.5 17))
(def v6 (double-array (range 6)))

(fvec/mult v6 25) ; scalar mult of a vector
(fvec/emult v6 v6) ; element-wise, i.e. Hadamard mult


(fmat/mulv (fmat/mat [[1 2]
                    [3 4]])
          (double-array [10 100])) ; a double-array vector is treated as column

;; This fails btw:
;(fmat/mulv (double-array [10 100])
;          (fmat/mat [[1 2]
;                    [3 4]]))


(comment
  ;; Doesn't work.
  ;; How do you make a column vector? Is it a column matrix?
(fmat/vtmul (double-array [10 100])
           (fmat/mat [[1 2]
                     [3 4]]))

;; No
(fmat/vtmul (fmat/mat [[10 100]])
           (fmat/mat [[1 2]
                     [3 4]]))

;; not this
(fmat/vtmul (Vec2. 10 100)
           (fmat/mat [[1 2]
                     [3 4]]))
)

;; OK so I don't know how to use vtmul, so the I can implement
;; fmat/vector mult either with vectors on right, or by using 
;; matrices to represent vectors, in which case I can use either
;; direction.
