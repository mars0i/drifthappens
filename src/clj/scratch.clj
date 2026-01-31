;; Initial ns statement copied from
;; https://generateme.github.io/fastmath/clay/vector_matrix.html#matrices
(ns scratch
  (:require ;[clojure.math :as m]
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

(comment
  (def m (fm/seq->double-double-array [[1.0M 2.0M][3.0M 4.0M]]))
  (type (fmat/entry m 1 1))
  (take-nth 5 (range 100))
)

(comment
  (def bindist (frand/distribution :binomial {:p 0.5 :trials 10}))
  (.probability bindist 5)
  (map (fn [n] (.probability bindist n)) (range 11))
)


(comment
  (defn umath/make-mat-powers
    [m expts]
    (let [size (first (fmat/shape m))]
      (loop [curr-mat m
             curr-expt 1
             remaining-expts expts
             acc-mats []]
        (let [next-mat (fmat/mulm m prev-mat)
              next-expt (inc curr-expt)
              rest-expts (rest expts)]
          (recur next-mat
                 (if (= curr-expt 
                        )))))))
)
