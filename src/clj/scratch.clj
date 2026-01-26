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


(defn plot-iterates
  "Plot lines connecting n iterates as y values with x representing number of
  iterations.  iterates should not be infinite."
  [iterates]
  (-> (tc/dataset {:x (range (count iterates))
                   :y iterates})
      ;(plotly/base {:=width 900})
      (plotly/layer-line {:=x :x, :=y, :y})
      plotly/plot
      (assoc-in [:data 0 :line :width] 1)))


;; pop size 50 with 50% A, 50% B.
(def init-freqs (wf/mkvec (concat (repeat 25 0.0) [1.0] (repeat 25 0.0))))

(def upmat (wf/right-mult-tran-mat 0.7 0.5 50))

(comment
  (fmat/mulv upmat init-freqs)
  (fmat/mulv (wf/mpow upmat 3) init-freqs)


)

(plot-iterates init-freqs)

(plot-iterates (fmat/mulv upmat init-freqs))
(plot-iterates (fmat/mulv (wf/mpow upmat 3) init-freqs))
(plot-iterates (fmat/mulv (wf/mpow upmat 8) init-freqs))
