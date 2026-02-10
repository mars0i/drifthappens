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
            [fastmath.random :as frand]
            ;[fastmath.dev.codox :as codox]
            ;[fastmath.dev.clay :as utls]
            ;[fastmath.dev.plotly :as plt]
            ;[fastmath.dev.ggplot :as gg]
            ;[fastmath.stats :as stats]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.tableplot.v1.plotly :as plotly]
            [tablecloth.api :as tc]
            [utils.fastmats :as mats]
            [utils.misc :as umisc])
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

;; ---

;; `sample-prob` computes
;; $p_i = \frac{iw_{A}}{w_{A}i + w_{B}(N-i)}\,,\,$
;; where $i=$`abs-freq-A`, $N=$`pop-size`, and $w__ =$`fit-_`.
(defn sample-prob
  "Returns the probability of sampling an individual of type A (rather than
  B) from a population of size pop-size, when exactly abs-freq-A members of
  the population have type A, and the fitnesses of A and B are fit-A and
  fit-B, respectively."
  [fit-A fit-B pop-size abs-freq-A]
  (let [abs-freq-B (- pop-size abs-freq-A)
        overall-fit-A (* abs-freq-A fit-A)
        overall-fit-B (* abs-freq-B fit-B)
        total-fit (+ overall-fit-A overall-fit-B)]
    (/ overall-fit-A total-fit)))

;; ---

;; TODO Why am I calculating this by hand rather than using frand/distribution?
;; e.g.:
;; ```
;; (def bindist (frand/distribution :binomial {:p 0.5 :trials 10}))
;; (map (fn [n] (.probability bindist n)) (range 11))
;; ```
;; `tran-prob` computes
;; $\binom{M}{j} p_i^j q_i^{M-j}$,
;; where $M=$`sample-size`, $j=$`As-sampled`, and $p_i$ is the value of values 
;; of `sample-prob` for $i=$`freq-A`, and $q_i=1-p_i$.
(defn tran-prob-obsolete
  "Returns the probability that sample of sample-size individuals will
  include exactly As-sampled A (not B) individuals, when the probability of
  choosing an A individual is sample-prob-A.  Note that sample-size can't
  be much more than 1000--else the function may return #NaN's or #Inf's.
  (big-tran-prob can handle larger sample sizes.)"
  [sample-prob-A sample-size As-sampled]
  (let [sample-prob-B (- 1 sample-prob-A)
        Bs-sampled (- sample-size As-sampled)]
    (* (fm/combinations sample-size As-sampled)
       (fm/pow sample-prob-A As-sampled)
       (fm/pow sample-prob-B Bs-sampled))))

;; TODO Why am I calculating this by hand rather than using frand/distribution?
;; e.g.:
;; ```
;; (def bindist (frand/distribution :binomial {:p 0.5 :trials 10}))
;; (map (fn [n] (.probability bindist n)) (range 11))
;; ```
;; `big-tran-prob` also computes $\binom{M}{j} p_i^j q_i^{M-j}$, but 
;; less precisely and it should allow $M\gg 1000$.  It uses the fact that
;$\binom{M}{j} p_i^j q_i^{M-j} = 
;;\exp \left[\ln\binom{M}{j} + \ln\left(p_i^j q_i^{M-j}\right)\right]$
;; However, in practice, if the sample size is much more than 1000, 
;; the values will all be zero.
(defn big-tran-prob-obsolete?
  "Returns the (approximate) probability that sample of sample-size
  individuals will include exactly As-sampled A (not B) individuals, when
  the probability of choosing an A individual is sample-prob-A.  This
  function can handle sample sizes much greater than 1000 (but the result
  may be that all values are zero)."
  [sample-prob-A sample-size As-sampled]
  (let [sample-prob-B (- 1 sample-prob-A)
        Bs-sampled (- sample-size As-sampled)
        log-combo (fm/log-combinations sample-size As-sampled)
        log-prob (fm/ln (* (fm/pow sample-prob-A As-sampled)
                           (fm/pow sample-prob-B Bs-sampled)))]
    (fm/exp (+ log-combo log-prob)))) ; replace `exp` with `qexp`?

;; ---

;; TODO NEED MORE TESTS OF THESE TWO VERSIONS

(defn tran-mat-elems
  "Creates a sequence of sequences of transition matrix elements that can
  be passsed to left-mult-tran-mat or right-mult-tran-mat.  The sum of
  the values in each inner sequence is (approximately) 1.0."
  [fit-A fit-B pop-size sample-size]
  (for [i (umisc/irange pop-size)]
    (let [sp (sample-prob fit-A fit-B pop-size i)
          bindist (frand/distribution :binomial {:p sp :trials sample-size})]
      (mapv (fn [n] (.probability bindist n)) (umisc/irange sample-size)))))

(defn tran-mat-elems-old
  "Creates a sequence of sequences of transition matrix elements that can
  be passsed to left-mult-tran-mat or right-mult-tran-mat.  The sum of
  the values in each inner sequence is (approximately) 1.0."
  [fit-A fit-B pop-size sample-size]
  (for [i (umisc/irange pop-size)]
    (let [sp (sample-prob fit-A fit-B pop-size i)]
      (for [j (umisc/irange sample-size)]
        (let [tp (tran-prob-obsolete sp sample-size j)]
          tp)))))

(comment
  (def dim 1000)
  (def m1 (mapv vec (tran-mat-elems 1.1 1.0 dim dim)))
  (def m2 (mapv vec (tran-mat-elems-old 1.1 1.0 dim dim)))
  (apply max (for [i (range dim), j (range dim)]
               (- ((m1 i) j) ((m2 i) j))))

  (def m1cat (apply concat m1))
  (apply min m1cat)
  (apply max m1cat)
  (double (/ (count (filter zero? m1cat))
             (count m1cat)))
)

(defn left-mult-tran-mat
  "Create a transition matrix in which the sum of values in each row is
  equal to 1.  Use this to multiply a row vector with the vector on the
  left and the matrix on the right."
  ([fit-A fit-B pop-size]
   (left-mult-tran-mat fit-A fit-B pop-size pop-size))
  ([fit-A fit-B pop-size sample-size]
   (mkmat (tran-mat-elems fit-A fit-B pop-size sample-size))))

(defn right-mult-tran-mat
  "Create a transition matrix in which the sum of values in each column is
  equal to 1.  Use this to multiply a column vector with the vector on the
  right and the matrix on the left."
  ([fit-A fit-B pop-size]
   (right-mult-tran-mat fit-A fit-B pop-size pop-size))
  ([fit-A fit-B pop-size sample-size]
   (mats/transpose
     (left-mult-tran-mat fit-A fit-B pop-size sample-size))))

;; ---
;; Tests and experiments for code above:

(comment
  (def tme (tran-mat-elems 1.0 1.0 1000 50))
  (map (partial apply +) tme) ; values should all be near 1.0
  (count tme) ; 1001
  (count (first tme)) ; 51
  (apply = (map count tme)) ; true
  (def lmtmraw (mkmat tme))
  (fmat/shape lmtmraw) ; [1001 51]
  (def lmtm (left-mult-tran-mat 1.0 1.0 1000 50))
  (fmat/shape lmtm) ; [1001 51]
  (def rmtm (right-mult-tran-mat 1.0 1.0 1000 50))
  (fmat/shape rmtm) ; [51 1001]
  (def mult-along-51 (fmat/mulm lmtm rmtm))
  (fmat/shape mult-along-51) ; [1001 1001]
  (def mult-along-1001 (fmat/mulm rmtm lmtm))
  (fmat/shape mult-along-1001)  ; [51 51]

  (def rmtm1K (right-mult-tran-mat 1.0 1.05 1000))
  (require '[utils.math :as um])
  (um/mpow rmtm1K 100)
)

(comment
  (defn compare-small-big
    [sample-prob-A sample-size As-sampled]
    (let [small (tran-prob sample-prob-A sample-size As-sampled)
          big (big-tran-prob sample-prob-A sample-size As-sampled)]
      (list small big (= small big) (- small big))))

  (compare-small-big 0.7 10 5)
  (compare-small-big 0.7 100 5)
  (compare-small-big 0.7 100 50)

  ;; These results are *really* small, often resolving to zero:

  (compare-small-big 0.7 1000 1)
  (compare-small-big 0.7 1000 5)
  (compare-small-big 0.7 1000 50)
  (compare-small-big 0.7 1000 350)
  (compare-small-big 0.7 1000 500)

  (compare-small-big 0.9 1000 5)
  (compare-small-big 0.9 1000 50)
  (compare-small-big 0.9 1000 350)
  (compare-small-big 0.9 1000 500)

  (compare-small-big 0.5 1000 5)
  (compare-small-big 0.5 1000 50)
  (compare-small-big 0.5 1000 350)
  (compare-small-big 0.5 1000 500)

  ;; These work in both versions:
  (compare-small-big 0.5 1050 5)
  (compare-small-big 0.5 1050 50)
  (compare-small-big 0.5 1050 350)
  ;; Here small-tran-prob returns #Inf,
  ;; while big-tran prob returns approximately 0.0075:
  (compare-small-big 0.5 1050 500)


  (compare-small-big 0.5 2000 5)
  (compare-small-big 0.5 2000 50)
  ;; small-tran-prob returns #NaN:
  (compare-small-big 0.5 2000 350)
  (compare-small-big 0.5 2000 500)

  (compare-small-big 0.5 10000 5)
  (compare-small-big 0.5 10000 50)
  ;; small-tran-prob returns #NaN, and big-tran-prob returns 0.0:
  (compare-small-big 0.5 10000 350)
  (compare-small-big 0.5 10000 500)
  (compare-small-big 0.5 10000 3500)
  (compare-small-big 0.5 10000 5000)

  ;; All of these values are zero, even in the middle.
  (def tran-probs-2K (map (partial big-tran-prob 0.5 2000) (range 2001)))
  (def tran-probs-2K (map (partial big-tran-prob 0.5 2000) (range 985 1015)))
  ;; So big-tran-prob isn't very useful.
  ;; Not sure it's useful to generate tran probs for N=1000, even with
  ;; small-tran-prob.  The probabilities are zero, but they're very small,
  ;; and sometimes zero:
  (def tran-probs-1K (map (partial small-tran-prob 0.5 1000) (range 1001)))
  (def tran-probs-1K (map (partial small-tran-prob 0.5 1000) (range 985 1025)))
  (def tran-probs-1K (map (partial small-tran-prob 0.5 1000) (range 5)))

  (fm/combinations 1000 500)
  (fm/exp (fm/log-combinations 1000 500))
  (fm/combinations 1000 5)
  (fm/combinations 1000 995)
  (fm/exp (fm/log-combinations 1000 5))
  (fm/exp (fm/log-combinations 1000 5))
  (fm/combinations 2000 1000)
  ;; Note that big-tran-prob works as well as it does because the
  ;; probability term reduces the size.  This is #Inf
  (fm/exp (fm/log-combinations 2000 1000))
  ;; because this value means that $e$ was raised to a high power:
  (fm/log-combinations 2000 1000)

  (small-tran-prob 0.5 500 250)
  (small-tran-prob 0.5 500 1)
  (small-tran-prob 0.5 500 499)
  (big-tran-prob 0.5 500 250)
  (big-tran-prob 0.5 500 1)
  (big-tran-prob 0.5 500 499)
)

