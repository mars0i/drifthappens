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

(def m33 (nc/mm 1.0 m32 m23))
(def m22 (nc/mm 1.0 m23 m32))

(def v (nn/dv [10 100 1000]))

(comment
  ;; Why are these returning (a copy of?) the original vector?
  ;; Well that's weird; now they're not. (Maybe I downloaded a bugfix?)
  (nc/mv m33 v) ; vector multiplication treats vector on right as column
  (nc/mv m23 v)

  (nc/mm m23 v) ; fails

  ;; weird naming
  (nc/mrows m23)
  (nc/ncols m23)
  (nc/mrows m32)
  (nc/ncols m32)
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


