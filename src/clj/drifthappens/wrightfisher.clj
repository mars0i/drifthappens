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
        Bs-sampled (- sample-size As-sampled)]
    (* (fm/combinations sample-size As-sampled)
       (fm/pow sample-prob-A As-sampled)
       (fm/pow sample-prob-B Bs-sampled))))

(defn tran-mat-elems
  [fit-A fit-B pop-size sample-size]
  (for [i (umisc/irange pop-size)]
    (let [sp (sample-prob fit-A fit-B pop-size i)]
      (for [j (umisc/irange sample-size)]
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

(comment
  (def tme (tran-mat-elems 1.0 1.0 1000 50))
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
)
)
