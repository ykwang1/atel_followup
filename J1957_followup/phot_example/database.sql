create database sentinel;
use sentinel;
create table watchlist(id varchar(32), ra float, decl float, field smallint, ccd tinyint, quad tinyint, filter varchar(3), last_updated float default 58787.5);
create table photdata (id varchar(32), mjd float, mag float, magerr float, fid tinyint, pid tinyint);
