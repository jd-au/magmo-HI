-- Table: public.field

-- DROP TABLE public.field;

CREATE TABLE public.field
(
    field_name character varying(30) COLLATE pg_catalog."default" NOT NULL,
    CONSTRAINT "field_pkey" PRIMARY KEY (field_name)
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.field
    OWNER to postgres;


-- Table: public.observation

-- DROP TABLE public.observation;

CREATE TABLE public.observation
(
    day integer NOT NULL,
    field_name character varying(30) COLLATE pg_catalog."default" NOT NULL,
    g_longitude double precision NOT NULL,
    g_latitude double precision NOT NULL,
    max_flux double precision,
    sn_ratio double precision,
    strong boolean,
    used boolean,
    CONSTRAINT observation_field_name FOREIGN KEY (field_name)
        REFERENCES public.field (field_name) MATCH SIMPLE
        ON UPDATE NO ACTION
        ON DELETE NO ACTION
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.observation
    OWNER to postgres;


-- Table: public.spectrum

-- DROP TABLE public.spectrum;

CREATE TABLE public.spectrum
(
    day integer NOT NULL,
    field_name character varying(30) COLLATE pg_catalog."default" NOT NULL,
    source character varying(100) COLLATE pg_catalog."default" NOT NULL,
    g_longitude double precision NOT NULL,
    g_latitude double precision NOT NULL,
    max_flux double precision,
    min_opacity double precision,
    max_opacity double precision,
    rms_opacity double precision,
    min_velocity double precision,
    max_velocity double precision,
    used boolean,
    opacity_range double precision,
    max_s_max_n double precision,
    continuum_sd double precision,
    max_em_std double precision,
    rating character varying(5) COLLATE pg_catalog."default",
    resolved boolean,
    filename character varying(1000) COLLATE pg_catalog."default",
    local_path character varying(1000) COLLATE pg_catalog."default",
    local_emission_path character varying(1000) COLLATE pg_catalog."default",
    CONSTRAINT spectrum_field_name FOREIGN KEY (field_name)
        REFERENCES public.field (field_name) MATCH SIMPLE
        ON UPDATE NO ACTION
        ON DELETE NO ACTION
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.spectrum
    OWNER to postgres;



-- Table: public.gas

-- DROP TABLE public.gas;

CREATE TABLE public.gas
(
    day integer NOT NULL,
    field_name character varying(30) COLLATE pg_catalog."default" NOT NULL,
    source character varying(100) COLLATE pg_catalog."default" NOT NULL,
    velocity double precision,
    em_velocity double precision,
    optical_depth double precision,
    temp_off double precision,
    temp_spin double precision,
    g_longitude double precision NOT NULL,
    g_latitude double precision NOT NULL,
    fwhm double precision,
    vel_diff double precision,
    equiv_width double precision,
    tau_v double precision,
    CONSTRAINT gas_field_name FOREIGN KEY (field_name)
        REFERENCES public.field (field_name) MATCH SIMPLE
        ON UPDATE NO ACTION
        ON DELETE NO ACTION
)
WITH (
    OIDS = FALSE
)
TABLESPACE pg_default;

ALTER TABLE public.gas
    OWNER to postgres;
