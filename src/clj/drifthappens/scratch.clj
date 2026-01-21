;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns drifthappens.scratch
  (:require ;[clojure.math :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [fastmath.core :as fm]
            ;[fastmath.dev.codox :as codox]
            ;[fastmath.dev.clay :as utls]
            ;[fastmath.dev.plotly :as plt]
            ;[fastmath.dev.ggplot :as gg]
            ;[fastmath.stats :as stats]
            [scicloj.kindly.v4.kind :as kind]
            [utils.misc :as misc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))


(def mat1 (mat/real-matrix [[1 2 3 4]
                            [5 6 7 8]]))

(def mat2 (mat/mat         [[1 2 3 4]
                            [5 6 7 8]]))

(def mat3 (mat/mat         [[1 2]
                            [5 6]]))
(type mat1)
(type mat2)
(type mat3)


(= mat1 mat2)

(mat/entry mat1 1 1)
(mat/entry mat2 1 1)

(mat/shape mat1)
(mat/shape mat2)

(def mm1 (mat/mulmt mat1 mat2))
(mat/shape mm1)
(mat/entry mm1 1 1)

(def mm2 (mat/mulm mat1 (mat/transpose mat2)))
(mat/shape mm2)
(mat/entry mm2 1 1)

(def mm3 (mat/mulm (mat/transpose mat1) mat2))
(mat/shape mm3)
(mat/entry mm3 1 1)
(mat/entry mm3 3 3)

(fm/combinations 4 2)


;; ---

(defn sample-prob
  [fit-A fit-B pop-size freq-A]
  (let [freq-B (- pop-size freq-A)
        overall-fit-A (* freq-A fit-A)
        overall-fit-B (* freq-B fit-B)
        total-fit (+ overall-fit-A overall-fit-B)]
    (/ overall-fit-A total-fit)))

(defn tran-prob
  [sample-prob-A sample-prob-B sample-size num-A]
  (let [num-B (- sample-size num-A)]
    (* (fm/combinations sample-size num-A)
       (fm/pow sample-prob-A num-A)
       (fm/pow sample-prob-B num-B))))

(defn tran-mat
  [fit-A fit-B pop-size sample-size]
  (let [pop-idxs (misc/irange pop-size)
        sample-idxs (misc/irange sample-size)
        sample-probs (map (partial sample-prob fit-A fit-B pop-size) pop-idxs) ; ???
        tran-probs (map identity sample-idxs)] ; FIXME
    ))

