(ns utils.math
  (:require [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]))

;; faster for larger m and n
;; uses mutual recursion, but number of calls is O(log2 n)
(defn mpow
  "Multiplies a square matrix by itself n times.  n must be a positive
  integer or a positive float with no fractional part."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond (not= h w) (do (print (str "mpow only multiplies square matrices; shape is [" h " " w "].")) nil) ; or throw?
          (or (<= n 0) (not (== (rem n 1) 0)))
          (do (print "mpow only accepts positive integer powers.") nil) ; or throw?
          :else (letfn [(either-pow [acc-mat k]
                          (cond (= k 1) acc-mat
                                (even? k) (even-pow acc-mat k)
                                :else (fmat/mulm m (even-pow acc-mat (dec k)))))
                        (even-pow [acc-mat k]
                          (let [acc-mat-part (either-pow acc-mat (/ k 2))]
                            (fmat/mulm acc-mat-part acc-mat-part)))]
                  (either-pow m (long n)))))) ; allow float integers

;; tail recursive, but slow for moderately large m and n
(defn mpow-slow
  "Multiplies a square matrix by itself n times."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond (not= h w) (do (print (str "mpow only multiplies square matrices; shape is [" h " " w "].")) nil) ; or throw?
          (or (<= n 0) (not (== (rem n 1) 0)))
          (do (print "mpow only accepts positive integer powers.") nil) ; or throw?
          :else (loop [i 1, acc-mat m]
                  (if (= i n)
                    acc-mat
                    (recur (inc i)
                           (fmat/mulm m acc-mat)))))))

(comment
  (def m (fm/seq->double-double-array [[1 2 3][4 5 6][7 8 9]]))
  (mpow m 10)
  (mpow-slow m 10)
  (mpow m 0)
  (mpow m 0.0)
  (mpow m -2)
  (mpow m 2.5)
  (mpow m 2.0) ; should succeed
  (def m1 (fm/seq->double-double-array [[1 2 3][4 5 6]]))
  (mpow m1 2)
)

