(ns neanderthaleg
  (:require [uncomplicate.neanderthal
             [core :as nc]
             [native :as nn]
             ;[real :as nreal]
             [random :as nrand]
            ]
            [uncomplicate.fluokitten.core :as ktn]
           ))



;; works
(def m23 (nn/dge 2 3 (range 6)))
(def m32 (nn/dge 3 2 (range 6)))

(def m33 (nc/mm 1.0 v32 v23))
(def m22 (nc/mm 1.0 v23 v32))

(def v (nn/dv [10 100 1000]))

(comment
  ;; Why are these returning (a copy of?) the original vector?
  (nc/mv m33 v)
  (nc/mv m23 v)

  (nc/mm m23 v) ; fails
)



(map identity (nn/dge 2 3 (range 6))) ;; lazy seq of 3 lazy pairs
(ktn/fmap identity (nn/dge 2 3 (range 6))) ; same 2x3 mat
(into [] (nn/dge 2 3 (range 6))) ; Clojure vector of lazy pairs

(comment
  ;; doesn't work?
  (nc/ge native-float 2 3 (range 6)) ;; worked once, but now it doesn't.  What was different?

  ;; don't work
  (ktn/fmap vector (nn/dge 2 3 (range 6))) ; error
  (ktn/fmap vec (nn/dge 2 3 (range 6))) ; error
  (into vec (nn/dge 2 3 (range 6))) ; error
)


