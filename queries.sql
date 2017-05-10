-- Multiple observations of the same field
select field_name, round(max(sn_ratio)*100)/100 as max_sn, round(min(sn_ratio)*100)/100, count(*)
from observation
where used = True
group by field_name
having count(*) > 1
order by 4 desc

-- Multiple spectra from the same field
select field_name, round(max(rms_opacity)*100)/100 as max_rms_opacity, round(min(rms_opacity)*100)/100 as min_rms_opacity, count(*)
from spectrum
group by field_name
order by 4 desc


-- List duplicate observations
select field_name, day, max_flux, round(sn_ratio*100)/100, strong, used from observation o
where field_name in (
select field_name
from observation
group by field_name
having count(*) > 1)
order by field_name, day