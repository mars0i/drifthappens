;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns scratch
  (:require ;[clojure.math :as m]
            [criterium.core :as crit]
            [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]
            [fastmath.random :as frand]
            [scicloj.kindly.v4.kind :as kind]
            [scicloj.tableplot.v1.plotly :as plotly]
            [tablecloth.api :as tc]
            [drifthappens.wrightfisher :as wf]
            [utils.plotly :as uplot]
            [utils.math :as umath]
            [utils.misc :as umisc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))


;; Experiment below shows with num-gens = 8 and half-interval = 8 that
;; choose-mat-powers-separately is faster than -sequentially (20 secs
;; rather than 26 secs).


(defn make-prob-states
  [tran-mats init-state]
  (mapv (fn [m] (fmat/mulv m init-state))
        tran-mats))

(def fit-B 1.0)
(def big-fit-A 1.01)  ; large size has sel
(def small-fit-A 1.0) ; small size is pure drift
(def interval 16) ; interval between generations to display
(def half-interval (/ interval 2))

(def increments (iterate (partial + 16) 0)) ; OLD VERSION
(def generations (map inc increments)) ; OLD VERSION
(def half-generations (map (fn [n] (inc (/ n 2))) increments)) ; OLD VERSION

(def num-gens 8)

;; These numbers need to be even:
(def small-N 10)
(def big-N 1000)

(def half-small-N (/ small-N 2))
(def half-big-N (/ big-N 2))

(def small-pop-init (wf/mkvec (concat (repeat half-small-N 0.0) [1.0] (repeat half-small-N 0.0))))
(comment
(def small-drift-mat (wf/right-mult-tran-mat small-fit-A fit-B (dec (count small-pop-init))))
(def small-tran-mats (doall (umath/choose-mat-powers-separately small-drift-mat (take num-gens generations)))) ; OLD VERSION
(def small-prob-states (make-prob-states small-tran-mats small-pop-init))
(def small-plots (mapv uplot/plot-both small-prob-states))
)

(def big-pop-init (wf/mkvec (concat (repeat half-big-N 0.0) [1.0] (repeat half-big-N 0.0))))
(comment
(def big-drift-mat (wf/right-mult-tran-mat big-fit-A fit-B (dec (count big-pop-init)))) ; use fit-B for fit-A to make them equal
(def big-tran-mats (doall (umath/choose-mat-powers-separately big-drift-mat (take num-gens generations)))) ; OLD VERSION
(def big-prob-states (make-prob-states big-tran-mats big-pop-init))
(def big-plots (mapv uplot/plot-lines big-prob-states))
)

;small-plots
;big-plots

(def N (dec (count big-pop-init)))
(def M (dec (count small-pop-init)))

;; Shrinking and expanding matrices.
;; Note M and N are swapped in the next two
(def predat-drift-mat (wf/right-mult-tran-mat small-fit-A fit-B N M)) ; 
(def reprod-drift-mat (wf/right-mult-tran-mat big-fit-A fit-B M N)) ; use fit-B for fit-A to make them equal
;; Combine into one square matrix:
(def pred-reprod-mat (fmat/mulm reprod-drift-mat predat-drift-mat))

;; We use half-generations here because each step involves two sampling processes.
;; So each generation is analogous to two generations in the small and big models.
(crit/bench
(def pred-reprod-tran-mats (doall (umath/choose-mat-powers-separately pred-reprod-mat (take num-gens half-generations)))) ; OLD VERSION
)

(comment
(def pred-reprod-prob-states (make-prob-states pred-reprod-tran-mats big-pop-init))
(def pred-reprod-plots (mapv uplot/plot-both pred-reprod-prob-states))

;pred-reprod-plots

;(def small-big-combo-plots (interleave small-plots big-plots pred-reprod-plots))

; RUN THE PLOTS:
;small-big-combo-plots 
)
