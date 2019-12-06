WITH
nearobj AS (
    SELECT
        nb.objID, r.RFGC
    FROM
        MyDB.RFGCfull AS r
        CROSS APPLY fGetNearbyObjEq(r.RAJ2000, r.DEJ2000, 0.16666666666666666) AS nb
    -- WHERE r.RFGC < 10
)
, panpos AS (
    SELECT
        o.objID, o.objName, o.raMean, o.decMean
    FROM
        MeanObjectView AS o
        INNER JOIN nearobj ON (o.objID=nearobj.objID)
)
, ser AS (
    SELECT s.*
    FROM
        StackModelFitSer AS s
        INNER JOIN nearobj ON (s.objID=nearobj.objID)
)
, seruniq AS (
    SELECT *
    FROM
    ( SELECT *
        ,ROW_NUMBER () OVER (
            PARTITION BY objID 
            ORDER BY bestDetection DESC, primaryDetection DESC
        ) AS rown -- enumerate inside common id group
        FROM ser
    ) AS t
    WHERE t.rown = 1 -- only first
)
, shape AS (
    SELECT s.*
    FROM ForcedGalaxyShape AS s
         INNER JOIN nearobj ON (s.objID=nearobj.objID)
)
, shapeuniq AS (
    SELECT *
    FROM ( 
        SELECT 
            *
          , ROW_NUMBER () OVER (
                PARTITION BY objID 
                ORDER BY
                CASE WHEN gGalFlags = 0 THEN 1
                     WHEN gGalFlags = 4 THEN 2
                     ELSE 3 END ASC,
                CASE WHEN rGalFlags = 0 THEN 1
                     WHEN rGalFlags = 4 THEN 2
                     ELSE 3 END ASC,
                CASE WHEN iGalFlags = 0 THEN 1
                     WHEN iGalFlags = 4 THEN 2
                     ELSE 3 END ASC,
                CASE WHEN zGalFlags = 0 THEN 1
                     WHEN zGalFlags = 4 THEN 2
                     ELSE 3 END ASC,
                CASE WHEN yGalFlags = 0 THEN 1
                     WHEN yGalFlags = 4 THEN 2
                     ELSE 3 END ASC
        ) AS rown -- enumerate inside common id group
        FROM shape
    ) AS t
    WHERE t.rown = 1 -- only first
)
, pet AS (
    SELECT p.*
    FROM
        StackPetrosian AS p
        INNER JOIN nearobj ON (p.objID=nearobj.objID)
)
, petuniq AS (
    SELECT *
    FROM
    ( SELECT *
        ,ROW_NUMBER () OVER (
            PARTITION BY objID 
            ORDER BY bestDetection DESC, primaryDetection DESC
        ) AS rown -- enumerate inside common id group
        FROM pet
    ) AS t
    WHERE t.rown = 1 -- only first
)
, kron AS (
    SELECT
        a.objID,
        a.bestDetection AS RadiusbestDetection, a.primaryDetection AS RadiusprimaryDetection,
        o.gMeanKronMag AS gkronMag, a.gKronRad AS gkronRadius,
        o.rMeanKronMag AS rkronMag, a.rKronRad AS rkronRadius,
        o.iMeanKronMag AS ikronMag, a.iKronRad AS ikronRadius,
        o.zMeanKronMag AS zkronMag, a.zKronRad AS zkronRadius,
        o.yMeanKronMag AS ykronMag, a.yKronRad AS ykronRadius
    FROM
        StackObjectAttributes AS a
        INNER JOIN nearobj ON (a.objID=nearobj.objID)
        LEFT JOIN MeanObjectView AS o on (o.objID=a.objID)
)
, kronuniq AS (
    SELECT *
    FROM
    ( SELECT *
        ,ROW_NUMBER () OVER (
            PARTITION BY objID 
            ORDER BY RadiusbestDetection DESC, RadiusprimaryDetection DESC
        ) AS rown -- enumerate inside common id group
        FROM kron
    ) AS t
    WHERE t.rown = 1 -- only first
)

SELECT 
    rfgc.*, 
    panpos.objID, 
    panpos.objName, 
    panpos.raMean, 
    panpos.decMean 
    ,
    seruniq.bestDetection AS serbestDetection,
    seruniq.primaryDetection AS serprimaryDetection
    ,
    seruniq.gSerRadius,seruniq.gSerAb,seruniq.gSerPhi,seruniq.gSerRa,seruniq.gSerDec,seruniq.gSerMag,
    seruniq.rSerRadius,seruniq.rSerAb,seruniq.rSerPhi,seruniq.rSerRa,seruniq.rSerDec,seruniq.rSerMag,
    seruniq.iSerRadius,seruniq.iSerAb,seruniq.iSerPhi,seruniq.iSerRa,seruniq.iSerDec,seruniq.iSerMag,
    seruniq.zSerRadius,seruniq.zSerAb,seruniq.zSerPhi,seruniq.zSerRa,seruniq.zSerDec,seruniq.zSerMag,
    seruniq.ySerRadius,seruniq.ySerAb,seruniq.ySerPhi,seruniq.ySerRa,seruniq.ySerDec,seruniq.ySerMag
    ,
    shapeuniq.gGalMinor,shapeuniq.gGalMajor,shapeuniq.gGalPhi,shapeuniq.gGalIndex,shapeuniq.gGalMag,shapeuniq.gGalFlags,
    shapeuniq.rGalMinor,shapeuniq.rGalMajor,shapeuniq.rGalPhi,shapeuniq.rGalIndex,shapeuniq.rGalMag,shapeuniq.rGalFlags,
    shapeuniq.iGalMinor,shapeuniq.iGalMajor,shapeuniq.iGalPhi,shapeuniq.iGalIndex,shapeuniq.iGalMag,shapeuniq.iGalFlags,
    shapeuniq.zGalMinor,shapeuniq.zGalMajor,shapeuniq.zGalPhi,shapeuniq.zGalIndex,shapeuniq.zGalMag,shapeuniq.zGalFlags,
    shapeuniq.yGalMinor,shapeuniq.yGalMajor,shapeuniq.yGalPhi,shapeuniq.yGalIndex,shapeuniq.yGalMag,shapeuniq.yGalFlags
    ,
    petuniq.bestDetection AS petbestDetection,
    petuniq.primaryDetection AS petprimaryDetection
    ,
    kronuniq.RadiusbestDetection AS kronRadiusbestDetection,
    kronuniq.RadiusprimaryDetection AS kronRadiusprimaryDetection
    , 
    petuniq.gpetMag,petuniq.gpetRadius,petuniq.rpetMag,petuniq.rpetRadius,petuniq.ipetMag,petuniq.ipetRadius,petuniq.zpetMag,petuniq.zpetRadius,petuniq.ypetMag,petuniq.ypetRadius, 
    kronuniq.gkronMag,kronuniq.gkronRadius,kronuniq.rkronMag,kronuniq.rkronRadius,kronuniq.ikronMag,kronuniq.ikronRadius,kronuniq.zkronMag,kronuniq.zkronRadius,kronuniq.ykronMag,kronuniq.ykronRadius 
INTO rfgc_nearby_multiband2

FROM MyDB.RFGCfull as rfgc
LEFT JOIN nearobj   ON nearobj.RFGC    = rfgc.RFGC
LEFT JOIN panpos    ON panpos.objID    = nearobj.objID
LEFT JOIN seruniq   ON seruniq.objID   = panpos.objID
LEFT JOIN shapeuniq ON shapeuniq.objID = seruniq.objID
LEFT JOIN kronuniq  ON kronuniq.objID  = shapeuniq.objID
LEFT JOIN petuniq   ON petuniq.objID   = kronuniq.objID

-- vim:ft=sqlanywhere
