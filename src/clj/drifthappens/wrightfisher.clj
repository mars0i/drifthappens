;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns drifthappens.wrightfisher
  (:require ;[clojure.math :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [fastmath.core :as fm]
            ;[fastmath.dev.codox :as codox]
            ;[fastmath.dev.clay :as utls]
            ;[fastmath.dev.plotly :as plt]
            ;[fastmath.dev.ggplot :as gg]
            ;[fastmath.stats :as stats]
            [scicloj.kindly.v4.kind :as kind]
            [tablecloth.api :as tc]
            [scicloj.tableplot.v1.plotly :as plotly]
            [utils.misc :as misc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))


(def mat1 (mat/real-matrix [[1 2 3 4]
                            [5 6 7 8]]))

(def mat2 (mat/mat         [[1 2 3 4]
                            [5 6 7 8]]))

(def mat3 (mat/mat         [[1 2]
                            [5 6]]))
(type mat1)
(type mat2)
(type mat3)


(= mat1 mat2)

(mat/entry mat1 1 1)
(mat/entry mat2 1 1)

(mat/shape mat1)
(mat/shape mat2)

(def mm1 (mat/mulmt mat1 mat2))
(mat/shape mm1)
(mat/entry mm1 1 1)

(def mm2 (mat/mulm mat1 (mat/transpose mat2)))
(mat/shape mm2)
(mat/entry mm2 1 1)

(def mm3 (mat/mulm (mat/transpose mat1) mat2))
(mat/shape mm3)
(mat/entry mm3 1 1)
(mat/entry mm3 3 3)

(fm/combinations 4 2)


;; ---

;; $$p_i = \frac{iw_{A}}{w_{A}i + w_{B}(N-i)}$$
;;
;; where $i=$`freq-A`, $N=$`pop-size`, and $w_\alpha =$`fit-α`.
(defn sample-prob
  [fit-A fit-B pop-size freq-A]
  (let [freq-B (- pop-size freq-A)
        overall-fit-A (* freq-A fit-A)
        overall-fit-B (* freq-B fit-B)
        total-fit (+ overall-fit-A overall-fit-B)]
    (/ overall-fit-A total-fit)))

;; $$\Pi_{ji} = \binom{M}{j} p_i^j q_i^{M-j}$$
;;
;; where $M=$`sample-size`, $j=$`num-A`, and $p_i$ and $q_i$ are
;; values of `sample-prob` for $i=$`freq-A`.
(defn tran-prob
  [sample-prob-A sample-size num-A]
  (let [sample-prob-B (- 1 sample-prob-A)
        num-B (- sample-size num-A)]
    (* (fm/combinations sample-size num-A)
       (fm/pow sample-prob-A num-A)
       (fm/pow sample-prob-B num-B))))

;; Better to multiply with a row vector on the left, since those 
;; vectors are more natural in Clojure and fastmath.
;; So (I think) it's row totals that sum to 1.  See Grimmett & Stirzacker
;; 3d p. 215.
;; Oh, but mulv mults by a col vec that's on the right. And vtmul
;; does a transpose first.  Need to learn aobut the vectors.

;; That means that the input is on the left, so num rows has to be
;; $N=$`pop-size`, while the output matches the column count,
;; $M=$`sample-size`.  So it's an $N\times M$ i.e.
;; `pop-size`$\times$`sample-size` matrix.
;;
;; i.e. for each $i\in 0, \ldots, N$ inclusive, generate a row with
;; tran prob values for each $j\in 0, \ldots, M$.

(defn tran-mat
  [fit-A fit-B pop-size sample-size]
  (let [pop-idxs (misc/irange pop-size)
        sample-idxs (misc/irange sample-size)
        rows (for [freq-A pop-idxs]
               (let [s-prob-A (sample-prob fit-A fit-B pop-size freq-A)]
                 (map (partial tran-prob s-prob-A sample-size)
                      sample-idxs)))]
    (mat/mat rows)))

(comment
  (def m (tran-mat 5 2 3 2))
  (mat/shape m)
;#object[org.apache.commons.math3.linear.Array2DRowRealMatrix 0x51bf848e
;        "Array2DRowRealMatrix
;        {{1.0,0.0,0.0},
;         {0.1975308642,0.4938271605,0.3086419753},
;         {0.0277777778,0.2777777778,0.6944444444},
;         {0.0,0.0,1.0}}"]

    (mat/mat->array2d m)

    ;; Should be a sequence of 1's:
    (map (fn [xs] (apply + xs))
         (mat/mat->array2d m))

)
