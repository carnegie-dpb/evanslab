--
-- count up the alt alleles per contig originating from a bunch of samples
--

DROP TABLE gbscounts;

CREATE TABLE gbscounts (
       contig          text     NOT NULL,
       loc             int      NOT NULL,
       ref             char(1)  NOT NULL,
       minmlog10p      double precision  NOT NULL DEFAULT 100.00,
       a               int      NOT NULL DEFAULT 0,
       c               int      NOT NULL DEFAULT 0,
       g               int      NOT NULL DEFAULT 0,
       t               int      NOT NULL DEFAULT 0
);

CREATE UNIQUE INDEX gbscounts_idx ON gbscounts(contig,loc);
