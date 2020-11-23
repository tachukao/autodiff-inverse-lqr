open Base
open Owl
module AD = Algodiff.D

let dir = Cmdargs.(get_string "-d" |> force ~usage:"-d [dir]")
let in_dir = Printf.sprintf "%s/%s" dir
let n = 2
let m = 2
let _T = 20
let n_samples = 20

(* true system *)
let a = Mat.of_arrays [| [| 1.; 1. |]; [| 0.; 1. |] |]
let b = Mat.eye m
let c = Mat.of_arrays [| [| 1.; 0. |] |]
let x0 = Mat.gaussian n n_samples
let q = Mat.(transpose c *@ c)
let r = Mat.of_arrays [| [| 0.1; 0.3 |] |] |> Mat.diagm

let k =
  let bt = Mat.transpose b in
  let p = Linalg.D.dare a b q r in
  let h = Mat.(bt *@ p *@ a) in
  let r_btb = Mat.(r + (bt *@ p *@ b)) in
  Linalg.D.linsolve r_btb h


let acl = Mat.(a - (b *@ k))

let run x0 =
  let rec run i xacc uacc x =
    if i < _T
    then (
      let u = Mat.(neg (k *@ x)) in
      let x = Mat.(acl *@ x) in
      run (Int.succ i) (x :: xacc) (u :: uacc) x)
    else
      ( xacc |> Array.of_list |> Arr.stack ~axis:0
      , uacc |> Array.of_list |> Arr.stack ~axis:0 )
  in
  run 0 [ x0 ] [] x0


(* true dynamics *)
let xs, us = run x0

(* now we optimise *)
module Prms = struct
  type 'a t = { q : 'a } [@@deriving prms]
end

open Prms
module O = Owl_opt_lbfgs.Make (Prms)

let unpack =
  let edn = AD.Mat.eye n in
  fun prms -> AD.Maths.((prms.q *@ transpose prms.q) + (F 1E-9 * edn))


let () =
  let a = AD.pack_arr a in
  let b = AD.pack_arr b in
  let r = AD.pack_arr r in
  let x0 = AD.pack_arr x0 in
  let xs = AD.pack_arr xs in
  let prms0 = { q = AD.Mat.eye n } in
  let f =
    let run acl k x0 =
      let rec run i xacc uacc x =
        if i < _T
        then (
          let u = AD.Maths.(neg k *@ x) in
          let x = AD.Maths.(acl *@ x) in
          run (Int.succ i) (x :: xacc) (u :: uacc) x)
        else
          ( xacc |> Array.of_list |> AD.Maths.stack ~axis:0
          , uacc |> Array.of_list |> AD.Maths.stack ~axis:0 )
      in
      run 0 [ x0 ] [] x0
    in
    fun _ prms ->
      let q = unpack prms in
      let k =
        let bt = AD.Maths.transpose b in
        let p = AD.Linalg.dare a b q r in
        let h = AD.Maths.(bt *@ p *@ a) in
        let r_btb = AD.Maths.(r + (bt *@ p *@ b)) in
        AD.Linalg.linsolve r_btb h
      in
      let acl = AD.Maths.(a - (b *@ k)) in
      let xs', _ = run acl k x0 in
      let xerr = AD.Maths.(l2norm_sqr' (xs' - xs)) in
      AD.Maths.(xerr / F (Float.of_int Int.(n_samples * _T)))
  in
  let s0 = O.init ~f ~prms0 () in
  let res = ref [] in
  O.min
    ~stop:(fun s ->
      let iter = O.iter s in
      let fv = O.fv s in
      let qhat = AD.unpack_arr (unpack O.(prms s)) in
      res := (fv, qhat) :: !res;
      Stdio.printf "%i %f\n%!" iter fv;
      false)
    s0
  |> ignore;
  let fvs, qhats = List.(rev !res) |> List.unzip in
  List.map ~f:(fun qhat -> Mat.(l2norm' (qhat - q)) |> Float.to_string) qhats
  |> Stdio.Out_channel.write_lines (in_dir "errors_in_q");
  fvs |> List.map ~f:Float.to_string |> Stdio.Out_channel.write_lines (in_dir "fvs")
