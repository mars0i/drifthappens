(ns utils.math
  (:require [clojure.math :as math]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]))

;; Most of matrix code has been moved elsewhere e.g. utils.fastmath

(def ln2 (math/log 2))

(defn log2
  "Returns the log base 2 of x."
  [x]
  (/ (math/log x) ln2))

(defn power-of-two?
  [x]
  (let [absl2 (log2 x)]
    (== absl2 (math/round absl2))))

(comment
  (power-of-two? 128)
)

