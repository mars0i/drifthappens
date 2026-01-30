(ns utils.math
  (:require [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]))

;; Fewer multiplications. Faster for larger m and n than simpler version.
;; Uses mutual recursion, but stack depth is O(log2 n).
(defn mpow
  "Multiplies a square matrix by itself n times.  n must be a positive
  integer or positive floating-point number with no fractional part."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond 
      (not= h w) 
      (print (str "mpow only multiplies square matrices. Matrix shape was [" h " " w "].")) ; throw?

      (or (<= n 0)
          (not (zero? (rem n 1))))
      (print (str "mpow only accepts positive integer powers. Exponent was " n)) ; throw?

      :else 
      (letfn [(either-pow [acc-mat k]
                (cond (= k 1) acc-mat
                      (even? k) (even-pow acc-mat k)
                      :else (fmat/mulm m (even-pow acc-mat (dec k)))))
              (even-pow [acc-mat k]
                (let [factor (either-pow acc-mat (/ k 2))]
                  (fmat/mulm factor factor)))]
        (either-pow m (long n)))))) ; allow float integers

;; Simple, tail recursive, but slow for moderately large m and n.
(defn mpow-slow
  "Multiplies a square matrix by itself n times.  n must be a positive
  integer or positive floating-point number with no fractional part."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond 
      (not= h w)
      (print (str "mpow only multiplies square matrices; shape is [" h " " w "].")) ; throw?

      (or (<= n 0) (not (== (rem n 1) 0)))
      (print "mpow only accepts positive integer powers.") ; throw?

      :else 
      (loop [i 1, acc-mat m]
        (if (= i n)
          acc-mat
          (recur (inc i)
                 (fmat/mulm m acc-mat)))))))

(defn make-mat-powers
  "Returns a sequence of transition matrices that are integer powers of of
  square matrix m for each exponent in expts.  Does not make use of smaller
  products previously produced."
  [m expts]
  (let [size (first (fmat/shape m))]
    (map (partial mpow m) expts)))

(defn mat-powers
  "Returns a lazy sequence of powers of a square matrixm, beginning with
  m^0, i.e. begining with the identity matrix of the same size."
  [m]
  (iterate (partial fmat/mulm m)
           (fmat/eye (first (fmat/shape m)))))

(comment
  (def m (fm/seq->double-double-array [[1 2 3][4 5 6][7 8 9]]))
  (mpow m 10)
  (mpow-slow m 10)
  (= (mpow m 10) (mpow-slow m 10))
  (mpow m 0)
  (mpow m 0.0)
  (mpow m -2)
  (mpow m 2.5)
  (mpow m 2.0) ; should succeed
  (def m1 (fm/seq->double-double-array [[1 2 3][4 5 6]]))
  (mpow m1 2)

  (mpow m 4)
  (mpow m 5)
  (make-mat-powers m [1 4 5])
  (def ms (mat-powers m))
  (take 8 m)

  (def m (fm/seq->double-double-array [[1 2 3][4 5 6][7 8 9]]))
  (def I (fmat/eye 3))
  (fmat/mulm I m)
  ; eval (effective-root-form): (fmat/mulm I m)
  ; (err) Execution error (ClassCastException) at fastmath.matrix.Mat3x3/mulm (matrix.clj:313).
  ; (err) class [[D cannot be cast to class fastmath.matrix.Mat3x3 ([[D is in module java.base of loader 'bootstrap'; fastmath.matrix.Mat3x3 is in unnamed module of loader clojure.lang.DynamicClassLoader @c1e1f4c)
  (def Ireal (fmat/eye 3 true))
  (fmat/mulm Ireal m)
  ; eval (effective-root-form): (fmat/mulm Ireal m)
  ; (err) Execution error (ClassCastException) at fastmath.matrix/eval41805$fn (matrix.clj:686).
  ; (err) class [[D cannot be cast to class org.apache.commons.math3.linear.RealMatrix ([[D is in module java.base of loader 'bootstrap'; org.apache.commons.math3.linear.RealMatrix is in unnamed module of loader 'app')
  ;; THIS WORKS:
  (def Idouble (fm/seq->double-double-array [[1 0 0][0 1 0][0 0 1]]))
  (fmat/mulm Idouble m)

  (fmat/diagonal 1 1 1)

  (type I)
  (type Ireal)
  (type Idouble)
)

