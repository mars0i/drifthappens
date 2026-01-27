Example illustrating an algorithm for matrix powers.

let mpow x r =
  let (frac_part, whole_part) = Pervasives.modf r in
  if frac_part <> 0. then failwith "mpow: fractional powers not implemented";
  let m, n = _matrix_shape x in assert (m = n);
  (* integer matrix powers using floats: *)
  if r < 1. then failwith "mpow: exponent is non-positive";
  let rec either_pow s acc =
     if s = 1. then acc
     else if mod_float s 2. = 0.
     then even_pow s acc
     else odd_pow s acc
  and even_pow s acc =
    let acc2 = either_pow (s /. 2.) acc in
    dot acc2 acc2
  and odd_pow s acc =
    dot x (even_pow (s -. 1.) acc)
  in either_pow r x
