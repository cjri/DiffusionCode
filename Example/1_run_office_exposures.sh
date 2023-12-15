for r in {1..500}
do
   echo "Generating exposure data for radius "$r
   ../run_diffusion --radius $r --opt_t0 0 --verb 0 --emm Cough --env Office
done
