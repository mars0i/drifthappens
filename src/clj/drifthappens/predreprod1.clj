;; Experiments that motivate claims in ~/papers/Drift20205/submissions/POBAM2026
;; about how the effects of e.g. an internal predation process that reduces N 
;; temporarily can lead to a great deal of drift even though N stays the same
;; from generation to generation.  Proof of concept.
;;
;; There are possible problems with generating the binomial
;; probabilities---they become small, are often zero--and with 
;; how long it takes to multiply the matrices when N is greater than
;; 1000 or 2000.  For the latter problem, consider using Neanderthal

(ns drifthappens.predreprod1
  (:require ;[clojure.math :as m]
            ;[criterium.core :as crit]
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


(defn make-prob-states
  [tran-mats init-state]
  (mapv (fn [m] (fmat/mulv m init-state))
        tran-mats))

(def fit-B 1.0)
(def big-fit-A 1.01)  ; large size has sel
(def small-fit-A 1.0) ; small size is pure drift
(def interval 16) ; interval between generations to display
(def half-interval (/ interval 2))

(def increments (iterate (partial + 16) 0))
(def generations increments)
(def half-generations (map (fn [n] (/ n 2)) increments))

(def num-gens 6)

;; These numbers need to be even:
(def small-N 10)
(def big-N 1000)

(def half-small-N (/ small-N 2))
(def half-big-N (/ big-N 2))


(def small-pop-init (wf/mkvec (concat (repeat half-small-N 0.0) [1.0] (repeat half-small-N 0.0))))
(def small-drift-mat (wf/right-mult-tran-mat small-fit-A fit-B (dec (count small-pop-init))))
;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def small-tran-mats (doall (umath/choose-mat-powers-separately small-drift-mat (take num-gens generations))))
(def small-prob-states (make-prob-states small-tran-mats small-pop-init))
(def small-plots (mapv uplot/plot-both small-prob-states))

;small-plots

(def big-pop-init (wf/mkvec (concat (repeat half-big-N 0.0) [1.0] (repeat half-big-N 0.0))))
(def big-drift-mat (wf/right-mult-tran-mat big-fit-A fit-B (dec (count big-pop-init)))) ; use fit-B for fit-A to make them equal
;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def big-tran-mats (doall (umath/choose-mat-powers-separately big-drift-mat (take num-gens generations))))
(def big-prob-states (make-prob-states big-tran-mats big-pop-init))
(def big-plots (mapv uplot/plot-lines big-prob-states))

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
;; NOTE: Consider replacing choose-mat-powers-separately with choose-mat-powers-sequentially if the exponents are closely spaced; it might be more efficient:
(def pred-reprod-tran-mats (doall (umath/choose-mat-powers-separately pred-reprod-mat (take num-gens half-generations))))
(def pred-reprod-prob-states (make-prob-states pred-reprod-tran-mats big-pop-init))
(def pred-reprod-plots (mapv uplot/plot-both pred-reprod-prob-states))

;pred-reprod-plots

(def small-big-combo-plots (interleave small-plots big-plots pred-reprod-plots))

; RUN THE PLOTS FROM ALL THREE POPULATIONS:
small-big-combo-plots 

