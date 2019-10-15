library( tidyverse )

synapser::synLogin()

syn <- synExtra::synDownloader( "~/data/DRIAD/dsRNA", ifcollision="overwrite.local" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

main <- function()
{
    vr <- c(`otssp-167` = "otssp167", `jq1` = "(+)-jq1" )
    
    ## Load the raw counts
    X <- openxlsx::read.xlsx( "190809_DGE_dsRNAmi_counts.xlsx" ) %>%
        as_tibble() %>% rename( Drug = Fluid.name, Count = Live.number.of.objects ) %>%
        mutate_at( "Drug", str_to_lower ) %>%
        mutate_at( "Drug", recode, !!!vr )

    ## Compute summary statistic and order accordingly
    v1 <- c( "lipofectamine", "empty", "dmso" )
    Xm <- X %>% group_by( Drug ) %>% summarize( Count = median(Count) ) %>%
        arrange( Count ) %>% mutate_at( "Drug", as_factor ) %>%
        mutate( Highlight = case_when(
                    Drug %in% v1 ~ "Control",
                    Drug == "dsrna_alone" ~ "dsRNA",
                    TRUE ~ "no") )

    ## Load the associated metadata
    M <- syn_csv( "syn11801537" ) %>% transmute( LINCSID = lincs_id, Drug = str_to_lower(name) )

    ## Load TAS vectors
    vt <- c("JAK1","JAK2","JAK3","TYK2", "SIK1", "SIK2", "SIK3")
    Y <- readRDS( "tas_vector_annotated_long.rds" ) %>%
        filter( fp_name == "morgan_normal" ) %>% select(data) %>%
        unnest() %>% filter( compound_id_source == "hmsl" ) %>%
        select( LINCSID=compound_id, Target=entrez_symbol, TAS=tas ) %>%
        filter( Target %in% vt ) %>% mutate_at( "Target", factor, vt )

    ## Combine everything into a single data frame
    vtas <- c("1 (<100 nM)" = "#b2182b",
              "2 (100-999 nM)" = "#ef8a62",
              "3 (1-10 uM)" = "#fddbc7",
              "10 (>10 uM)" = "grey85")
    Z <- X %>% select(Drug) %>% distinct %>%
        left_join(M, by="Drug") %>% inner_join(Y, by="LINCSID") %>%
        mutate_at( "TAS", recode, !!!set_names(names(vtas), c(1,2,3,10)) ) %>%
        mutate_at( "TAS", as_factor ) %>%
        mutate_at( "Drug", factor, levels(Xm$Drug) )

    ## Plotting elements
    pal <- c(Control="steelblue", dsRNA="tomato", no="gray")
    gl <- guide_legend(nrow=2, byrow=TRUE, reverse=TRUE)
    etxt <- partial( element_text, face="bold" )
    theme_bold <- function() {theme( axis.title=etxt(size=11), axis.text=etxt(size=9) )}
    
    ## Individual points and summary statistic
    g1 <- ggplot( Xm, aes(x=Drug, y=Count) ) + theme_bw() + coord_flip() +
        geom_bar( aes(fill=Highlight), stat="identity", alpha=0.75 ) + 
        scale_fill_manual( values=pal, guide=FALSE ) +
        geom_point( data=X ) + theme_bold() + ylab("Nuclei Count") +
        scale_y_continuous(position="right")

    ## Overall distribution of values
    g2 <- ggplot( Xm, aes(x=Count) ) + theme_bw() +
        xlab( "Nuclei Count" ) + xlim( layer_scales(g1)$y$range$range ) +
        geom_density( adjust=.3 ) + ylab( "Density" ) +
        geom_vline( data=filter(Xm, Highlight != "no"),
                   aes(xintercept=Count, color=Highlight) ) +
        scale_color_manual( values=pal, guide=FALSE ) + theme_bold()

    ## Tas scores
    g3 <- ggplot( Z, aes(y=Drug, x=Target, fill=TAS) ) + theme_bw() +
        geom_tile() + scale_y_discrete(drop=FALSE, position="right") +
        scale_x_discrete(position="top") +
        scale_fill_manual( values = vtas, guide=gl ) + theme_bold() +
        theme( legend.position="bottom", axis.title.y=element_blank(),
              axis.text.x=element_text(angle=90, vjust=.5),
              legend.title=etxt(size=12), legend.text=etxt(size=10) )
    g3l <- cowplot::get_legend(g3)
    g3 <- g3 + guides( fill=FALSE )

    gg <- egg::ggarrange( g1, g3, g2, ncol=2, heights=c(9,1), widths=c(3,1) )
    ggp <- cowplot::ggdraw(gg) + cowplot::draw_plot( g3l, 0.33, -0.42 )
    ggsave( ggp, filename="dsRNA.pdf", width=9.5, height=11 )
}

