;; ### Wright-Fisher models
;; drifthappens/wrightfisher.clj

;; See tips/fastmathMatices.clj for tips on using fastmath matrices and vectors.

;; Note in "left/right multiplication" of a matrix and a vector, "right" or 
;; "left" refers to the location of the vector.

;; ---

;; (ns drifthappens.wrightfisher ...)
^:kindly/hide-code
;; Initial ns statement copied from https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns drifthappens.wrightfisher
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
            [scicloj.tableplot.v1.plotly :as plotly]
            [tablecloth.api :as tc]
            [utils.misc :as um])
  (:import [fastmath.vector ArrayVec Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))

;; ---


;; generateme/fastmath supports a variety of different vector and
;; matrix representations, but they don't all have the same conveniences
;; and don't necessarily work with the multiplication operators.
;; These are aliases for the current best choice for me.  See
;; https://clojurians.zulipchat.com/#narrow/channel/151924-data-science/topic/Matrix-vector.20multiplication.20in.20fastmath/near/569973078
(def mkvec double-array)
(def mkmat fm/seq->double-double-array)

;; TODO ? This could be made more efficient for large powers by
;; computing half the power and then multiplying the result
;; by itself. cf. owl_linalg_generic.ml in the OCaml Owl lib.
;; TODO Move somewhere else.
(defn mpow
  "Multiply a square matrix by itself n times."
  [m n]
  (if (or (<= n 0) (not (== (rem n 1) 0)))
    (do (print "mpow only accepts positive integer powers.")
        nil) ; or throw?
    (loop [i 1, acc-mat m]
      (if (= i n)
        acc-mat
        (recur (inc i)
               (fmat/mulm m acc-mat))))))

(comment

  (fvec/pow [1 2 3] 2)

  (== 0 0.0 0N 0M 0/1)

  (def m42 (fmat/mat [[2 3]
                      [4 5]]))

  (fmat/mulm m42 m42) ; {{16.0,21.0},{28.0,37.0}}"]
  (fmat/mulm m42 (fmat/mulm m42 m42)) ; {{116.0,153.0},{204.0,269.0}}"]
  (fmat/mulm m42 
             (fmat/mulm m42
                        (fmat/mulm m42 m42))) ; {{844.0,1113.0},{1484.0,1957.0}}"]
  (mpow m42 0)
  (mpow m42 -1)
  (mpow m42 1.00002)
  (mpow m42 1)
  (mpow m42 2)
  (mpow m42 3)
  (mpow m42 4)

  (def unit (fmat/mat [[2]]))
  (fmat/mat->array2d unit) ; row vec
)

;; ---

;; $p_i = \frac{iw_{A}}{w_{A}i + w_{B}(N-i)}\,,\,$
;; where $i=$`abs-freq-A`, $N=$`pop-size`, and $w__ =$`fit-_`.
(defn sample-prob
  "Returns the probability of sampling an individual of type A (rather than
  B) from a population of size pop-size, when exactly abs-freq-A members of
  the population have type A, and the fitnesses of A and B are fit-A and
  fit-B, respectively."
  [fit-A fit-B pop-size abs-freq-A]
  (let [num-B (- pop-size abs-freq-A)
        overall-fit-A (* abs-freq-A fit-A)
        overall-fit-B (* num-B fit-B)
        total-fit (+ overall-fit-A overall-fit-B)]
    (/ overall-fit-A total-fit)))

;; $\Pi_{ji} = \binom{M}{j} p_i^j q_i^{M-j}\,,\,$
;; where $M=$`sample-size`, $j=$`As-sampled`, and $p_i$ is the value of values 
;; of `sample-prob` for $i=$`freq-A`, and $q_i=1-p_i$.
;; 
;; NOTE: fm/combinations, and hence the following function, returns #NaN if 
;; $N=$`sample-size` gets much above 1000.
;; This seems to be due to the binomial coefficient getting really small.
;; So maybe this isn't Clojure's or fastmath's problem.  You'd need to use
;; some kind of approximation in any language, I think.  (Fastmath has a
;; log-combinations.  Could that be helpful?)
(defn tran-prob
  "Returns the probability that sample of sample-size individuals will
  include exactly As-sampled A (not B) individuals, when the probability of
  choosing an A individual is sample-prob-A.  Note that sample-size can't
  be much more than 1000--else the function will return #NaN."
  [sample-prob-A sample-size As-sampled]
  (let [sample-prob-B (- 1 sample-prob-A)
        num-B (- sample-size As-sampled)]
    (* (fm/combinations sample-size As-sampled)
       (fm/pow sample-prob-A As-sampled)
       (fm/pow sample-prob-B num-B))))

(comment
  (defn tran-probs
    "Create a sequence of transition probabilities for each number of A's
    between 0 and sample-size.  (Summing them should be very close to 1.)"
    [sample-prob-A sample-size]
    (map (partial tran-prob sample-prob-A sample-size)
         (um/irange sample-size)))

  (tran-probs 0.6 4.0)
  (reduce + (tran-probs 0.6 4.0))
  (reduce + (tran-probs 0.5 1000))
  (fm/combinations 1025 500)
  (tran-probs 0.4 200)
)

(defn tran-mat-elems
  [fit-A fit-B pop-size sample-size]
  (for [i (um/irange pop-size)]
    (let [sp (sample-prob fit-A fit-B pop-size i)]
      (for [j (um/irange sample-size)]
        (let [tp (tran-prob sp sample-size j)]
          tp)))))

(defn left-mult-tran-mat
  "Create a transition matrix in which the sum of values in each row is
  equal to 1.  Use this to multiply a row vector with the vector on the
  left and the matrix on the right."
  ([fit-A fit-B pop-size]
   (left-mult-tran-mat fit-A fit-B pop-size pop-size))
  ([fit-A fit-B pop-size sample-size]
   (mkmat (tran-mat-elems fit-A fit-B pop-size sample-size)))) ; each row sum = 1

(defn right-mult-tran-mat
  "Create a transition matrix in which the sum of values in each column is
  equal to 1.  Use this to multiply a column vector with the vector on the
  right and the matrix on the left."
  ([fit-A fit-B pop-size]
   (right-mult-tran-mat fit-A fit-B pop-size pop-size))
  ([fit-A fit-B pop-size sample-size]
   (fmat/transpose
     (left-mult-tran-mat fit-A fit-B pop-size sample-size))))

