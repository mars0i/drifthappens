(ns utils.math
  (:require [fastmath.vector :as fvec]
            [fastmath.matrix :as fmat]
            [fastmath.core :as fm]))

;; Fewer multiplications. Faster for larger m and n than simpler version.
;; Uses mutual recursion, but stack depth is O(log2 n).
(defn mpow
  "Multiplies a square matrix by itself n times.  n must be a nonnegative
  integer or positive floating-point number with no fractional part. Returns
  the identity nmatrix of the same shape when n=0."
  [m n]
  (let [[h w] (fmat/shape m)]
    (cond 
      (not= h w) 
      (print (str "mpow only multiplies square matrices. Matrix shape was [" h " " w "].")) ; throw?

      (or (< n 0) (not (zero? (rem n 1))))
      (print (str "mpow only accepts nonnegative integer powers. Exponent was " n)) ; throw?

      (zero? n)
      (fmat/mat->array2d (fmat/eye h)) ; at present need to convert eye to Java Array2D for mult with Array2D

      :else 
      (letfn [(either-pow [acc-mat k]
                (cond (= k 1) acc-mat
                      (even? k) (even-pow acc-mat k)
                      :else (fmat/mulm m (even-pow acc-mat (dec k)))))
              (even-pow [acc-mat k]
                (let [factor (either-pow acc-mat (/ k 2))]
                  (fmat/mulm factor factor)))]
        (either-pow m (long n)))))) ; allow float integers

(defn mat-powers
  "Returns a lazy sequence of powers of a square matrix beginning with
  m^0, i.e. begining with the identity matrix of the same size.  (Specific
  elements can be extracted e.g. using take-nth, keep-indexed, or 
  utils/misc/elems-at.)"
  [m]
  (iterate (partial fmat/mulm m)
           (fmat/mat->array2d ; at present need to convert eye to Java Array2D for mult with Array2D
             (fmat/eye (first (fmat/shape m))))))

(defn choose-mat-powers-sequentially
  "Returns a lazy sequence of transition matrices that are integer powers
  of square matrix m for each exponent in expts (which may be in any
  order). Exponents must be nonnegative integers, which may be in a float
  representation."
  [m expts]
  (let [expts-set (set expts)
        n (count expts-set)] ; don't count dupes
    (take n (keep-indexed (fn [i m] (when (expts-set i) m))
                          (mat-powers m)))))

;; FIXME floats are causing this to not return
(defn choose-mat-powers-separately
  "Returns a lazy sequence of transition matrices that are integer powers
  of square matrix m for each exponent in expts (which may be in any
  order). Exponents must be nonnegative integers, which may be in a float
  representation."
  [m expts]
  (map (partial mpow m) (map long expts)))

(comment
  (def m (fm/seq->double-double-array [[1 2 3][4 5 6][7 8 9]]))
  (def ms (mat-powers m))
  (take 10 ms)
  (nth ms 9)
  (nth ms 10)
  (def mpsep (choose-mat-powers-separately m [0 1 4 9]))
  (def mpsep2 (choose-mat-powers-separately m [0.0 1.0 4.0 9.0]))
  (def mpseq (choose-mat-powers-sequentially m [0 1 4 9]))
  (def mpseq2 (choose-mat-powers-sequentially m [1.0 4.0 9])) ; loops forever why?
)


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


(comment
;; There's a bug here that causes an infinite loop though the second
;; clause of the second if.  I don't understand it, but another verison
;; using `keep-indexed` works, so I'm not going to bother to fix this
;; unless there's a good reason to resuscitate it.
(defn make-mat-powers-sequentially-bad
  [m expts]
  (let [powseq (mat-powers m)]
    (loop [i 0, pows (mat-powers m), es expts, acc-pows []]
      (if (empty? expts)
        acc-pows
        (if (= i (first es))
          (recur (inc i)
                 (rest es)
                 (rest pows)
                (conj acc-pows (first pows)))
          (recur (inc i)
                 es ; don't move on until we find one
                 (rest pows)
                 acc-pows))))))
)

(comment
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
)

