;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns scratch
  (:require ;[clojure.math :as m]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.tableplot.v1.plotly :as plotly]
            [tablecloth.api :as tc]
            [drifthappens.wrightfisher :as wf]
            [utils.misc :as misc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))


(defn plot-lines
  "Plot lines connecting values in ys as y values with x representing the
  index of the element in ys.  Lines are drawn from point to point."
  [ys]
  (-> (tc/dataset {:x (range (count ys))
                   :y ys})
        (plotly/layer-line  {:=x :x, :=y, :y});)
      plotly/plot
      (assoc-in [:data 0 :line :width] 1)))

(defn plot-dots
  "Plot values in ys with dots t as y values with x representing the
  index of the element in ys."
  [ys]
  (-> (tc/dataset {:x (range (count ys))
                   :y ys})
        (plotly/layer-point {:=x :x, :=y, :y});)
      plotly/plot
      (assoc-in [:data 0 :line :width] 1)))

(defn plot-both
  "Plot lines connecting values in ys as y values with x representing the
  index of the element in ys.  Dots are plotted for poijnts, and lines are
  drawn from point to point."
  [ys]
  (-> (tc/dataset {:x (range (count ys))
                   :y ys})
        (plotly/layer-line  {:=x :x, :=y, :y});)
        (plotly/layer-point {:=x :x, :=y, :y});)
      plotly/plot
      (assoc-in [:data 0 :line :width] 1)))

;(comment
  (let [init-freqs (wf/mkvec (concat (repeat 500 0.0) [1.0] (repeat 500 0.0)))
        upmat (wf/right-mult-tran-mat 0.99 1.0 (dec (count init-freqs)))]
    (plot-lines init-freqs)
    (plot-lines (fmat/mulv upmat init-freqs))
    (plot-lines (fmat/mulv (wf/mpow upmat 2) init-freqs))
    ;(plot-lines (fmat/mulv (wf/mpow upmat 3) init-freqs))
    (plot-lines (fmat/mulv (wf/mpow upmat 8) init-freqs))
    (plot-lines (fmat/mulv (wf/mpow upmat 16) init-freqs))
    ;(plot-lines (fmat/mulv (wf/mpow upmat 40) init-freqs))
    (plot-lines (fmat/mulv (wf/mpow upmat 100) init-freqs))
  )
;)

(comment
  (let [init-freqs (wf/mkvec (concat (repeat 50 0.0) [1.0] (repeat 50 0.0)))
        upmat (wf/right-mult-tran-mat 1.0 1.0 (dec (count init-freqs)))]
    (kind/fragment
      [
       (plot-dots init-freqs)
       (plot-dots (fmat/mulv upmat init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 2) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 3) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 8) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 16) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 30) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 40) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 50) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 60) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 70) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 80) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 90) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 100) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 200) init-freqs))
       (plot-dots (fmat/mulv (wf/mpow upmat 500) init-freqs))
      ]
    )
  )
)
