(ns neanderthaleg
  (:require [uncomplicate.neanderthal
             [core :as nc]
             [native :as nn]
             ;[real :as nreal]
             ;[random :as nrand]
            ]
            [uncomplicate.fluokitten.core :as ktn]
           ))



(def v (nn/dv (range 50)))
(v 49)

;; works
(nn/dge 2 3 (range 6))
(nn/dge 3 2 (range 6))
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


