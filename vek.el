;;; vek.el --- opinionated geometric vector library. -*- lexical-binding: t; -*-

;; Author: Ryan Burnside
;; Version: 1.0.0
;; Keywords: math geometry vector vectors

;;; Commentary:

;; This provides geometric vector operations

;;; Code:

;; Do not use vector literals, modifying them is undefined!
;; Prefer vek

(defun vek (n &rest nums)
  "Returns a vector. Ensures all values are numberp."
  (let* ((vals (cons n nums)))
    (assert (every #'numberp vals))
    (apply #'vector vals)))

;; Macros enable /direct/ access by place
(defmacro vek-x (v)
  "setf-able synonym for aref v 0."
  `(aref ,v 0))

(defmacro vek-y (v)
  "setf-able synonym for aref v 1."
  `(aref ,v 1))

(defmacro vek-z (v)
  "setf-able synonym for aref v 2."
  `(aref ,v 2))

(defun veknp (v length)
  "Predicate to ensure a vector is both a vector and has length v"
  (and (vectorp v)
       (every #'numberp v)
       (= (length v) length)))

(defun vek2p (v)
  "Predicate to see if a vector is a 2 element vector."
  (veknp v 2))

(defun vek3p (v)
  "Predicate to see if a vector is a 3 element vector."
  (veknp v 3))

(defun vek4p (v)
  "Predicate to see if a vector is a 4 element vector."
  (veknp v 4))

(defun vek-zeros (size)
  "Return a vector of 0.0s with length of size."
  (make-vector size 0.0))

(defun vek-zeros-p (v)
  "Predicate for a vector full of 0."
  (and (vectorp v)
       (every #'zerop v)))

(defun vek-ones (size)
  "return a vector of 1.0s with a length of size."
  (make-vector size 1.0))

(defun vek-ones-p (v)
  "Predicate for a vector full of 1."
  (and (vectorp v)
       (every (apply-partially #'= 1) v)))

(defun vek-coerce (v size)
  "Promote or demote a vector to a size."
  (let ((delta (- size (length v))))
    (cond ((zerop delta) v)
          ((plusp delta) (vconcat v (vek-zeros (abs delta))))
          ((minusp delta) (subseq v 0 size)))))

(defun vek-coerce-all (v &rest vecs)
  "Promote all vectors to the dimension of the longest vector."
  (let* ((all-vecs (cons v vecs))
         (max-length (apply #'max
                            (mapcar #'length all-vecs))))
    (mapcar (lambda (vek)
              (vek-coerce vek max-length))
            all-vecs)))

(defun vek-map-elts (fn &rest vecs)
  "Perform a map operation on vectors promoting the shorter ones."
  (apply #'map 'vector fn (apply #'vek-coerce-all vecs)))

(defun vek+ (v &rest vecs)
  "Summation of atleast 1 vector."
  (apply #'vek-map-elts #'+ (cons v vecs)))

(defun vek- (v &rest vecs)
  "Subtraction of atleast 1 vector."
  (apply #'vek-map-elts #'- (cons v vecs)))

(defun vek* (v scale)
  "Multiply a vector by a scaler."
  (map 'vector
       (apply-partially #'* scale)
       (lambda (n) (* n scale))
       v))

(defun vek/ (v scale)
  "Divide a vector by a scaler."
  (map 'vector
       (lambda (n) (/ n (float scale)))
       v))

(defun vek-dot (v v2) ;; TODO make &rest version
  "Find the dot product of vectors."
  (cl-loop for n across v
           for n2 across v2
           sum (* n n2)))

(defun vek-cross (v v2)
  "Find the cross product of two vectors."
  (assert (and (vek3p v)
               (vek3p v2)))
  (vek (- (* (vek-y v) (vek-z v2)) (* (vek-z v) (vek-y v2)))
       (- (* (vek-z v) (vek-x v2)) (* (vek-x v) (vek-z v2)))
       (- (* (vek-x v) (vek-y v2)) (* (vek-y v) (vek-x v2)))))

(defun vek-negate (v)
  "Multiplies all components by -1."
  (vek- v))

(defun vek-mag (v)
  "Find magnitude of a vector."
  (sqrt (apply #'+ (mapcar (lambda (n)
                             (* n n))
                           v))))

(defun vek-gauss (mu sigma)
  "Generate a random number from a Gaussian distribution with mean MU and standard deviation SIGMA."
  (let* ((u1 (cl-random 1.0))
         (u2 (cl-random 1.0))
         (theta (* 2.0 pi (cl-random 1.0)))
         (z (sqrt (- 1.0 (* 2.0 (log u1)))))
         (x (* z (cos theta)))
         (y (* z (sin theta))))
    (+ mu (* sigma x))))

(defun vek-random (dims)
  "Return a random unit vector.
   We don't use rand for each element because that bunches at the corners."
  (assert (plusp dims))
  (let* ((vec (cl-loop repeat dims
                       collect (vek-gauss 0 1)))
         (mag (float (expt (apply #'+
                                  (mapcar (lambda (n)
                                            (* n n))
                                          vec))
                           .5))))
    (if (zerop mag)
        (let ((v (vek-zeros dims)))
          (setf (vek-x v) 1.0)
          v)
      (apply #'vector (cl-loop for i in vec
                               collect (/ i mag))))))

(defun vek-angle2D (v)
  (atan (vek-y v) (vek-x v)))

;; TODO
;; (vek-stretch v magnitude)
;; (vek-rotate2D v delta)

(provide 'vek)
;;; vek.el ends here
