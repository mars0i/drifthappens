(ns utils.misc)

;; By John Collins at https://stackoverflow.com/a/68476365/1455243
(defn irange
  "Inclusive range function: end element is included."
  ([start end step]
   (take-while (if (pos? step) #(<= % end) #(>= % end)) (iterate #(+ % step) start)))
  ([start end]
   (irange start end 1))
  ([end]
   (irange 0 end))
  ([] (range)))


;; Surprisingly, Clojure doesn't seem to have this built in, though
;; select-keys is close.  Consider using take-nth if you want every nth
;; element.
(defn elems-at
  "Return a sequence containing elements of xs at indexes in idxs,
  which must be sorted in increasing order."
  [xs idxs]
  (loop [acc [], ys xs, is idxs, i 0]
    (cond
      (or (empty? ys) (empty? is)) acc
      (= i (first is)) (recur (conj acc (first ys))
                              (rest ys)
                              (rest is) ; drop curr index only if used
                              (inc i))
      :else (recur acc, (rest ys), is, (inc i)))))

(comment
  (def xs (range 100 200 10))
  (elems-at xs [1 4 3])
  (elements-at xs [1 3 4 9])
)

;; This version is inefficient if xs is not an array.
;(defn elems-at
;  [xs idxs]
;  (map (partial nth xs) idxs))
