;; ### Wright-Fisher models
;; drifthappens.wrightfisher, defined in drifthappens/wrightfisher.clj

^:kindly/hide-code
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

;; See tips/fastmathMatices.clj for tips on using fastmath matrices
;; and vectors.


;; ---

;; $p_i = \frac{iw_{A}}{w_{A}i + w_{B}(N-i)}\,,\,$
;; where $i=$`freq-A`, $N=$`pop-size`, and $w__ =$`fit-_`.
(defn sample-prob
  "Returns the probability of sampling an individual of type A (rather than
  B) from a population of size pop-size, when exactly num-A members of the
  population have type A, and the fitnesses of A and B are fit-A and fit-B,
  respetively."
  [fit-A fit-B pop-size num-A]
  (let [num-B (- pop-size num-A)
        overall-fit-A (* num-A fit-A)
        overall-fit-B (* num-B fit-B)
        total-fit (+ overall-fit-A overall-fit-B)]
    (/ overall-fit-A total-fit)))

;; $\Pi_{ji} = \binom{M}{j} p_i^j q_i^{M-j}\,,\,$
;; where $M=$`sample-size`, $j=$`num-A`, and $p_i$ is the value of values 
;; of `sample-prob` for $i=$`freq-A`, and $q_i=1-p_i$.
(defn tran-prob
  "Returns the probability that sample of sample-size individuals will
  include exactly num-A A (not B) individuals, when the probability of
  choosing an A individual is sample-prob-A."
  [sample-prob-A sample-size num-A]
  (let [sample-prob-B (- 1 sample-prob-A)
        num-B (- sample-size num-A)]
    (* (fm/combinations sample-size num-A)
       (fm/pow sample-prob-A num-A)
       (fm/pow sample-prob-B num-B))))

(defn tran-probs
  "Create one row/column of transition probabilities that sum to 1."
  [sample-prob-A sample-size num-A]
  )

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
