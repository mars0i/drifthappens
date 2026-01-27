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
            [utils.plotly :as up]
            [utils.math :as umath]
            [utils.misc :as umisc])
  (:import [fastmath.vector Vec2 Vec3 Vec4]
           [fastmath.matrix Mat2x2 Mat3x3 Mat4x4]))


;; ??: 
;;  <groupId>org.apache.commons</groupId>
;;     <artifactId>commons-math3</artifactId>
;;     <version>3.6.1</version>
;; </dependency>
;; 
;; Once set up, we can use the RealMatrix interface and its
;; Array2DRowRealMatrix implementation to create our usual matrices.
;; The constructor of the implementation class takes a
;; two-dimensional double array as its parameter:
;; 
;; RealMatrix matrix = new Array2DRowRealMatrix(/* a two dimensions double array */);

;; https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/linear/Array2DRowRealMatrix.html
;; https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/linear/AbstractRealMatrix.html


(def init-freqs (wf/mkvec (concat (repeat 500 0.0) [1.0] (repeat 500 0.0))))
(def upmat (wf/right-mult-tran-mat 0.95 1.0 (dec (count init-freqs))))
;(type upmat)

(up/plot-lines init-freqs)
(up/plot-lines (fmat/mulv (umath/mpow upmat 2) init-freqs))
(up/plot-lines (fmat/mulv (umath/mpow upmat 4) init-freqs))
(up/plot-lines (fmat/mulv (umath/mpow upmat 20) init-freqs))
(up/plot-lines (fmat/mulv (umath/mpow upmat 40) init-freqs))
;(up/plot-lines (fmat/mulv (umath/mpow upmat 80) init-freqs))

(comment
  (def upmat2 (umath/mpow upmat 2))
  (def upmat4 (umath/mpow upmat2 2))
  (def upmat16 (umath/mpow upmat4 2))
  (def upmat256 (umath/mpow upmat16 2))
  (def upmat65K (umath/mpow upmat256 2))
  (up/plot-lines init-freqs)
  (up/plot-lines (fmat/mulv upmat init-freqs))
  (up/plot-lines (fmat/mulv upmat2 init-freqs))
  (up/plot-lines (fmat/mulv upmat4 init-freqs))
  (up/plot-lines (fmat/mulv upmat16 init-freqs))
  (up/plot-lines (fmat/mulv upmat256 init-freqs))
  (up/plot-lines (fmat/mulv upmat65K init-freqs))
)

(comment
  (let [init-freqs (wf/mkvec (concat (repeat 50 0.0) [1.0] (repeat 50 0.0)))
        upmat (wf/right-mult-tran-mat 1.0 1.0 (dec (count init-freqs)))]
    (kind/fragment
      [
       ;(up/plot-dots init-freqs)
       ;(up/plot-dots (fmat/mulv upmat init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 2) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 3) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 8) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 16) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 30) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 40) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 50) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 60) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 70) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 80) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 90) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 100) init-freqs))
       (up/plot-dots (fmat/mulv (umath/mpow upmat 200) init-freqs))
       ;(up/plot-dots (fmat/mulv (umath/mpow upmat 500) init-freqs))
      ]
    )
  )
)
