(ns drifthappens.predreprod1
  (:require ;[clojure.math :as m]
            [criterium.core :as crit]
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

(def increments (iterate (partial + 16) 0)) ; OLD VERSION
(def generations (map inc increments)) ; OLD VERSION
(def half-generations (map (fn [n] (inc (/ n 2))) increments)) ; OLD VERSION

(def num-gens 15)

;; These numbers need to be even:
(def small-N 10)
(def big-N 500)

(def half-small-N (/ small-N 2))
(def half-big-N (/ big-N 2))

(def small-pop-init (wf/mkvec (concat (repeat half-small-N 0.0) [1.0] (repeat half-small-N 0.0))))
(def small-drift-mat (wf/right-mult-tran-mat small-fit-A fit-B (dec (count small-pop-init))))

;; I think the make-mat-powers version is faster than the iterate version, even though it has to
;; redo matrix multiplications, because it saves multiplications by aggregating mults and then 
;; reusing them within a call to mpow.  I think this is because I'm only
;; taking every 16th matrix. If I was using every one of them, then mat-powers would probably be faster.
;;
;; (What would really be cool is some way of combining them so that you could reuse the calcs that mpow
;; does, across calls to it.  (I've haven't found memoize to be useful in the past.))
(def small-tran-mats (doall (take num-gens (take-nth interval (umath/mat-powers small-drift-mat))))) ; NEW VERSION: 20% slower (!)
;(def small-tran-mats (doall (umath/make-mat-powers small-drift-mat (take num-gens generations)))) ; OLD VERSION

(comment ; for testing
(def small-prob-states (make-prob-states small-tran-mats small-pop-init))
(def small-plots (mapv uplot/plot-both small-prob-states))
) ; comment for testing

(def big-pop-init (wf/mkvec (concat (repeat half-big-N 0.0) [1.0] (repeat half-big-N 0.0))))
(def big-drift-mat (wf/right-mult-tran-mat big-fit-A fit-B (dec (count big-pop-init)))) ; use fit-B for fit-A to make them equal
(crit/bench
(def big-tran-mats (doall (take num-gens (take-nth interval (umath/mat-powers big-drift-mat))))) ; NEW VERSION
)
(crit/bench
(def big-tran-mats (doall (umath/make-mat-powers big-drift-mat (take num-gens generations)))) ; OLD VERSION
)
(comment ; for testing
(def big-prob-states (make-prob-states big-tran-mats big-pop-init))
(def big-plots (mapv uplot/plot-lines big-prob-states))
) ; comment for testing

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

;; We use half-generations here because each step involves two sampling
;; processes.  So each generation is analogous to two generations in the
;; small and big models.
(def pred-reprod-tran-mats (doall (take num-gens (take-nth half-interval (umath/mat-powers pred-reprod-mat))))) ; NEW VERSION
;(def pred-reprod-tran-mats (doall (umath/make-mat-powers pred-reprod-mat (take num-gens half-generations)))) ; OLD VERSION
(comment ; for testing
(def pred-reprod-prob-states (make-prob-states pred-reprod-tran-mats big-pop-init))
(def pred-reprod-plots (mapv uplot/plot-both pred-reprod-prob-states))

;pred-reprod-plots

(def small-big-combo-plots (interleave small-plots big-plots pred-reprod-plots))

; RUN THE PLOTS:
small-big-combo-plots 
) ; comment for testing


(comment
  (fmat/shape predat-drift-mat)
  (fmat/shape reprod-drift-mat)
  (fmat/shape pred-reprod-mat)
  (map fmat/shape pred-reprod-tran-mats)
)
