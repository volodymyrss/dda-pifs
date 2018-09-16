import subprocess
import dataanalysis as da
import dataanalysis 
import ddosa
import glob
import ast
import pyfits,pywcs
import os

class ii_pif(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_cat=ddosa.CatExtract
    #input_image=ddosa.ii_skyimage
    input_ic=ddosa.IBIS_ICRoot
    input_bins=ddosa.ImageBins
    input_gti=ddosa.ibis_gti

    copy_cached_input=False

    cached=True

    #version="thermfov"
    #version="std"
    version="oof100"

   # write_caches=[dataanalysis.TransientCache,ddosa.MemCacheIntegralFallback,ddosa.MemCacheIntegralIRODS]
   # read_caches=[dataanalysis.TransientCache,ddosa.MemCacheIntegralFallback,ddosa.MemCacheIntegralFallbackOldPath,ddosa.MemCacheIntegralIRODS]

    off_edge_pixels=None

    #ii_pif_binary=os.environ['COMMON_INTEGRAL_SOFTDIR']+"/ii_pif/ii_pif_oof/ii_pif"
    ii_pif_binary=os.environ['COMMON_INTEGRAL_SOFTDIR']+"/i_one/ii_pif_therm_outfov/ii_pif"
    #ii_pif_binary="ii_pif"

    def get_cat(self):
        return self.input_cat.cat.get_path()
    
    def main(self):
        try:
            catpath=self.get_cat()
        except Exception as e:
            print "an exception getting cat:",e
            self.empty_results=True
            return

        ddosa.construct_gnrl_scwg_grp(self.input_scw,[\
                    #self.input_cat.cat.get_path(),
                    self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]",
                    self.input_gti.output_gti.get_path()
                ])

        ddosa.import_attr(self.input_scw.scwpath+"/swg.fits",["OBTSTART","OBTEND","TSTART","TSTOP","SW_TYPE","TELAPSE","RA_SCX","DEC_SCX","RA_SCZ","DEC_SCZ"])
        ddosa.set_attr({'ISDCLEVL':"BIN_I"})
        ddosa.set_attr({'INSTRUME':"IBIS"},"og.fits")

        ddosa.construct_gnrl_scwg_grp_idx([\
                    "og.fits",
                ])
        ddosa.set_attr({'ISDCLEVL':"BIN_I"},"og_idx.fits")

        ddosa.construct_og([\
                    "og_idx.fits",
                ])
        ddosa.set_attr({'ISDCLEVL':"BIN_I"},"ogg.fits")


        ddosa.remove_withtemplate("isgri_model.fits(ISGR-PIF.-SHD.tpl)")

        ht=ddosa.heatool(self.ii_pif_binary)
        #ht=ddosa.heatool("ii_pif")
        ht['inOG']=""
        ht['outOG']="ogg.fits[1]"
        ht['inCat']=catpath
        ht['mask']=self.input_ic.ibisicroot+"/mod/isgr_mask_mod_0003.fits[ISGR-MASK-MOD,1,IMAGE]"
#        ht['deco']=self.input_ic.ibisicroot+"/mod/isgr_deco_mod_0008.fits[ISGR-DECO-MOD,1,IMAGE]"
        ht['tungAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0010.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['aluAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0011.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['leadAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0012.fits[ISGR-ATTN-MOD,1,BINTABLE]"
 #       ht['covrMod']=self.input_ic.ibisicroot+"/mod/isgr_covr_mod_0002.fits[1]"
        ht['num_band'] = len(self.input_bins.bins)
        ht['E_band_min'] = " ".join([str(a[0]) for a in self.input_bins.bins])
        ht['E_band_max'] = " ".join([str(a[1]) for a in self.input_bins.bins])
        if self.off_edge_pixels is not None:       
            ht['AllowOffEdge'] = self.off_edge_pixels
        ht.run()

        self.pifs=da.DataFile("isgri_model.fits")
       # self.skyima=DataFile("isgri_sky_ima.fits")
       # self.skyres=DataFile("isgri_sky_res.fits")

class CatFromImaging(ddosa.DataAnalysis):
    input_cat=ddosa.ii_skyimage

    def main(self):
        skyres=pyfits.open(self.input_cat.skyres.get_path())[2].data # eband?

        fn="isgri_cat_from_image.fits"
        tpl="ISGR-SRCL-CAT.tpl"
        ddosa.remove_withtemplate(fn+"("+tpl+")")

        ht=ddosa.heatool("dal_create")
        ht['template']=tpl
        ht['obj_name']=fn
        ht.run()

        fo=pyfits.open(fn)

        fo['ISGR-SRCL-CAT'] = pyfits.BinTableHDU.from_columns(fo['ISGR-SRCL-CAT'].columns, nrows=len(skyres),header=fo['ISGR-SRCL-CAT'].header)
        fo['ISGR-SRCL-CAT'].data['RA_OBJ']=skyres['RA_FIN']
        fo['ISGR-SRCL-CAT'].data['DEC_OBJ']=skyres['DEC_FIN']
        fo.writeto(fn,clobber=True)

        self.cat=da.DataFile(fn)


class ii_pif_fromimaging(ii_pif):
    input_cat=ddosa.ii_skyimage
    #input_cat=CatFromImaging

    def get_cat(self):
        return self.input_cat.srclres.get_path()      


class evts_extract(ddosa.DataAnalysis):
    input_events=ddosa.ISGRIEvents
    input_scw=ddosa.ScWData
    input_cat=ddosa.CatExtract
    #input_image=ddosa.ii_skyimage
   # input_ic=ddosa.IBIS_ICRoot
  #  input_bins=ddosa.ImageBins
    input_pifs=ii_pif
    input_gti=ddosa.ibis_gti
    input_dead=ddosa.ibis_dead

    cached=True

    def main(self):
        att=self.input_scw.auxadppath+"/attitude_historic.fits"
        if os.path.exists(att) or os.path.exists(att+".gz"):
            att=self.input_scw.auxadppath+"/attitude_historic.fits[AUXL-ATTI-HIS,1,BINTABLE]"
            attp=att
        else:
            att=self.input_scw.auxadppath+"/attitude_snapshot.fits[AUXL-ATTI-SNA,1,BINTABLE]"
            attp_g=glob.glob(self.input_scw.auxadppath+"/attitude_predicted_*.fits*")
            attp=attp_g[0]+"[AUXL-ATTI-PRE,1,BINTABLE]"
        
        orb=self.input_scw.auxadppath+"/orbit_historic.fits"
        if os.path.exists(orb) or os.path.exists(orb+".gz"):
            orb=self.input_scw.auxadppath+"/orbit_historic.fits[AUXL-ORBI-HIS,1,BINTABLE]"
            orbp=orb
        else:
            orb=self.input_scw.auxadppath+"/orbit_snapshot.fits[AUXL-ORBI-SNA,1,BINTABLE]"
            orbp_g=glob.glob(self.input_scw.auxadppath+"/orbit_predicted_*.fits*")
            orbp=orbp_g[0]+"[AUXL-ORBI-PRE,1,BINTABLE]"

        ddosa.construct_gnrl_scwg_grp(self.input_scw,[\
                    #self.input_cat.cat.get_path(),
                    self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]",
                    self.input_scw.scwpath+"/isgri_events.fits[ISGR-EVTS-ALL]", \
                   # self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]", \
                    self.input_gti.output_gti.get_path(),
                    self.input_events.events.get_path(),
                    attp,
                    orbp,
                    self.input_dead.output_dead.get_path(),
                    self.input_pifs.pifs.get_path()
                ])


        ddosa.import_attr(self.input_scw.scwpath+"/swg.fits",["OBTSTART","OBTEND","TSTART","TSTOP","SW_TYPE","TELAPSE","RA_SCX","DEC_SCX","RA_SCZ","DEC_SCZ"])
        ddosa.set_attr({'ISDCLEVL':"BIN_I"})
        ddosa.set_attr({'INSTRUME':"IBIS"},"og.fits")

        ddosa.construct_gnrl_scwg_grp_idx([\
                    "og.fits",
                ])
        ddosa.set_attr({'ISDCLEVL':"BIN_I"},"og_idx.fits")

        ddosa.construct_og([\
                    "og_idx.fits",
                ])
        ddosa.set_attr({'ISDCLEVL':"BIN_I"},"ogg.fits")

        source_evts="source_evts.fits"
        ddosa.remove_withtemplate(source_evts+"(ISGR-EVTS-ALL.tpl)")

        ee=ddosa.heatool("evts_extract")
        ee['group']="ogg.fits" 
        ee['events']=source_evts
        ee['instrument']="IBIS" 
        ee['sources']=self.input_cat.cat.get_path()
        ee['gtiname']="MERGED_ISGRI" 
        ee['pif']='yes'
        ee['pifDOL']=self.input_pifs.pifs.get_path()
        ee['deadc']='yes'
        ee['attach']="no"
        ee['barycenter']=1
        ee['timeformat']=0
        ee['instmod']=""
        ee.run()

