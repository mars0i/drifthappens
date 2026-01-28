(ns utils.math
  (:require [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]))

;; Fewer multiplications. Faster for larger m and n than simpler version.
;; Uses mutual recursion, but stack depth is O(log2 n).
(defn mpow
  "Multiplies a square matrix by itself n times.  n must be a positive
  integer or positive floating-point number with no fractional part."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond 
      (not= h w) 
      (print (str "mpow only multiplies square matrices. Matrix shape was [" h " " w "].")) ; throw?

      (or (<= n 0)
          (not (zero? (rem n 1))))
      (print (str "mpow only accepts positive integer powers. Exponent was " n)) ; throw?

      :else 
      (letfn [(either-pow [acc-mat k]
                (cond (= k 1) acc-mat
                      (even? k) (even-pow acc-mat k)
                      :else (fmat/mulm m (even-pow acc-mat (dec k)))))
              (even-pow [acc-mat k]
                (let [factor (either-pow acc-mat (/ k 2))]
                  (fmat/mulm factor factor)))]
        (either-pow m (long n)))))) ; allow float integers

;; Simple, tail recursive, but slow for moderately large m and n.
(defn mpow-slow
  "Multiplies a square matrix by itself n times.  n must be a positive
  integer or positive floating-point number with no fractional part."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond 
      (not= h w)
      (print (str "mpow only multiplies square matrices; shape is [" h " " w "].")) ; throw?

      (or (<= n 0) (not (== (rem n 1) 0)))
      (print "mpow only accepts positive integer powers.") ; throw?

      :else 
      (loop [i 1, acc-mat m]
        (if (= i n)
          acc-mat
          (recur (inc i)
                 (fmat/mulm m acc-mat)))))))

(comment
  (def m (fm/seq->double-double-array [[1 2 3][4 5 6][7 8 9]]))
  (mpow m 10)
  (mpow-slow m 10)
  (= (mpow m 10) (mpow-slow m 10))
  (mpow m 0)
  (mpow m 0.0)
  (mpow m -2)
  (mpow m 2.5)
  (mpow m 2.0) ; should succeed
  (def m1 (fm/seq->double-double-array [[1 2 3][4 5 6]]))
  (mpow m1 2)
)

