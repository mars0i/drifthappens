(ns utils.math
  (:require [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]))

;; tail recursive, but slow for large m and n
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

;; Uses mutual recursion, but number of calls is O(log2 n).
(defn mpow
  "Multiplies a square matrix by itself n times."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond (not= h w) (do (print (str "mpow only multiplies square matrices; shape is [" h " " w "].")) nil) ; or throw?
          (or (<= n 0) (not (== (rem n 1) 0)))
          (do (print "mpow only accepts positive integer powers.") nil) ; or throw?
          :else (letfn [(either-pow [acc-mat k]
                          (cond (= k 1) acc-mat
                                (even? k) (even-pow acc-mat k)
                                :else (odd-pow acc-mat k)))
                        (even-pow [acc-mat k]
                          (let [acc-mat-part (either-pow acc-mat (/ k 2))]
                            (fmat/mulm acc-mat-part acc-mat-part)))
                        (odd-pow [acc-mat k]
                          (fmat/mulm m (even-pow acc-mat (dec k))))]
                  (either-pow m n)))))

(comment
  (def m1 (fm/seq->double-double-array [[1 2 3][4 5 6]]))
  (fmat/shape m1)
  (def m (fm/seq->double-double-array [[1 2 3][4 5 6][7 8 9]]))
  (time (mpow m 100))
  (time (mpow-slow m 100))
)

