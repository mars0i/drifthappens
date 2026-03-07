(ns scratch
  (:require ;[clojure.math :as m]
            ;[criterium.core :as crit]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]
            [fastmath.random :as frand]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.tableplot.v1.plotly :as plotly]
            [tablecloth.api :as tc]
            [drifthappens.wrightfisher :as wf]
            [utils.plotly :as uplot]
            [utils.fastmats :as mats]
            [utils.misc :as umisc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))


(defn binomdist
  [N p]
  (frand/distribution :binomial {:p p :trials N}))

(def maxN 300)

(def N (double 10))

(def mybinom (binomdist N 0.5))

(def yo
  (concat
    (map (fn [x] (.probability mybinom x)) (umisc/irange N))
    (repeat (- maxN N) 0)))


(uplot/plot-dots-lines yo)


