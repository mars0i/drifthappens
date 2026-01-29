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
            [utils.plotly :as uplot]
            [utils.math :as umath]
            [utils.misc :as umisc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))

(defn make-tran-mats
  [m expts]
  (let [size (first (fmat/shape m))]
    (doall
      (cons (fmat/eye size)
            (map (partial umath/mpow m)
                 expts)))))

(defn make-prob-states
  [tran-mats init-state]
  (mapv (fn [m] (fmat/mulv m init-state))
        tran-mats))


;; pop size 50 with 50% A, 50% B.
(def small-pop-init (wf/mkvec (concat (repeat 25 0.0) [1.0] (repeat 25 0.0))))
(def small-drift-mat (wf/right-mult-tran-mat 1.0 1.0 (dec (count small-pop-init))))
(def small-tran-mats (make-tran-mats small-drift-mat [1 5 10 15 20 25 30 35 40 45 50 55 60]))
(def small-prob-states (make-prob-states small-tran-mats small-pop-init))
(def small-plots (mapv uplot/plot-both small-prob-states))


;; pop size 500 with 50% A, 50% B.
(def big-pop-init (wf/mkvec (concat (repeat 250 0.0) [1.0] (repeat 250 0.0))))
(def big-drift-mat (wf/right-mult-tran-mat 1.0 1.0 (dec (count big-pop-init))))
(def big-tran-mats (make-tran-mats big-drift-mat [1 5 10 15 20 25 30 35 45 50 55 60]))
(def big-prob-states (make-prob-states big-tran-mats big-pop-init))
(def big-plots (mapv uplot/plot-lines big-prob-states))

;small-plots
;big-plots

(def N (dec (count big-pop-init)))
(def M (dec (count small-pop-init)))

(def predat-drift-mat (wf/right-mult-tran-mat 1.0 1.0 N M))
(def reprod-drift-mat (wf/right-mult-tran-mat 1.0 1.0 M N))
(def pred-reprod-mat (fmat/mulm reprod-drift-mat predat-drift-mat))

(def pred-reprod-tran-mats (make-tran-mats pred-reprod-mat [1 5 10 15 20 25 30 35 45 50 55 60]))
(def pred-reprod-prob-states (make-prob-states pred-reprod-tran-mats big-pop-init))
(def pred-reprod-plots (mapv uplot/plot-both pred-reprod-prob-states))

;pred-reprod-plots

(def small-big-combo-plots (interleave small-plots big-plots pred-reprod-plots))

small-big-combo-plots 


(comment
  (fmat/shape predat-drift-mat)
  (fmat/shape reprod-drift-mat)
  (fmat/shape pred-reprod-mat)
  (map fmat/shape pred-reprod-tran-mats)
)
