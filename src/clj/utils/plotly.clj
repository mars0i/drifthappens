^:kindly/hide-code
(ns utils.plotly
  (:require [clojure.math :as m] ; new in Clojure 1.11
            [scicloj.tableplot.v1.plotly :as plotly]
            [scicloj.kindly.v4.kind :as kind]
            ;[scicloj.kindly.v4.api :as kindly]
            [tablecloth.api :as tc]
            [utils.misc :as um]
            ))


(def plot-steps 400)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; CONVENIENCE PLOTTING FUNCTIONS

(defn plot-lines
  "Plot lines connecting values in ys as y values with x representing the
  index of the element in ys.  Lines are drawn from point to point."
  [ys]
  (kind/fragment
    [(kind/md "Yow!")
     (-> (tc/dataset {:x (range (count ys))
                      :y ys})
         (plotly/layer-line  {:=x :x, :=y, :y})
         plotly/plot
         (assoc-in [:data 0 :line :width] 1))]))

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn assoc-into-trace
  "Given a Plotly plot EDN, and an index for one of its traces (i.e. one of
  the maps in the vector that's the value of :data), Uses assoc-in to add
  or replace the value of the sequence of keys in ks with v."
  [plot trace-idx ks v]
  (assoc-in plot
            (into [:data trace-idx] ks)
            v))

(comment
  (def fake-plot {:data [{:x [1 2 3],
                          :y [4 5 6],
                          :xaxis3 {:anchor "y1"}}
                         {:x [10 20 30],
                          :y [40 50 60],
                          :xaxis {:anchor "y1"}}]
                  :layout {}})

  (-> fake-plot
      (assoc-into-trace 0 [:xaxis3 :anchor] "x4")
      (assoc-into-trace 0 [:xaxis3 :domain] [0.3 0.95]))
)

(defn assoc-into-traces
  "Applies the same insertion/addition to several traces using
  assoc-into-trace. trace-idxs is a sequence of indexes into the traces in
  the vector that's the value of :data in plot. See that assoc-into-trace
  for more info."
  [plot trace-idxs ks v]
  (reduce (fn [p idx] (assoc-into-trace p idx ks v))
          plot trace-idxs))

(defn set-subplot-order
  [plot trace-idxs subplot-order]
  (let [x-order-str (str "x" subplot-order)
        y-order-str (str "y" subplot-order)]
    (-> plot
        (assoc-into-traces trace-idxs [:xaxis] x-order-str)
        (assoc-into-traces trace-idxs [:yaxis] y-order-str))))

;; -----------------------------------------------

;; [xmin x-max] is the first arg, for use with -> , because it's likely
;; that different fns will share same x coords.
(defn fn2dataset
  "Generates a TMD dataset of points in which the coordinates are (x, f(x))
  from x-min to x-max, with an additional field named catkey, with a single
  value catval.  This key/val pair can be used to identify these coordinates
  after they are merged into another dataset."
  [[x-min x-max] catkey catval f & {:keys [steps]}] 
  (let [x-increment (/ (double (- x-max x-min))
                       (or steps plot-steps))
        xs (um/irange x-min x-max x-increment)
        ys (map f xs)]
    (-> {:x xs, :y ys, catkey catval}
        tc/dataset)))

(defn equalize-display-units
  "Given a Tableplot plot (in either Hanami-key form or full Plotly EDN),
  adds Plotly settings to force the displayed x and y units to be equal.
  (Note if the plot's display dimensions are too large in one dimension,
  there will be extra space in the plot outside of the plotted data.)"
  [plot]
  (-> plot
      plotly/plot
      (assoc-in [:layout :yaxis :scaleanchor] :x)
      (assoc-in [:layout :yaxis :scaleratio] 1)))

(defn set-line-width
  "Sets the ith Plotly plot in :data to have width (default: 2).  (Note
  that you can set the default line width by including the key `:=mark-sze`
  with a numeric value in the map passed to plotly/layer-line.)"
  [plot i width]
  (assoc-in plot [:data i :line :width] width))

(defn set-line-dash
  "Sets the ith Plotly plot in :data to have dash value dashed which should
  be a string. See
  https://plotly.com/javascript/reference/scatter/#scatter-line-dash for
  possible values."
  [plot i dashed]
  (assoc-in plot [:data i :line :dash] dashed))
